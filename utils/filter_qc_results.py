#!/usr/bin/env python
import argparse
import json
import os


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )

    p.add_argument("--in", dest="in_dir", required=True,
                   help="Input directory with all subjects/bundles.")
    p.add_argument("--qc", dest="qc_report", required=True,
                   help="JSON QC report.")
    p.add_argument("--out", required=True,
                   help="Output directory")
    p.add_argument("--force_remove", nargs="+",
                   help="List of bundles to force remove even when the qc "
                        "report says to keep it.")

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if not os.path.exists(args.in_dir):
        raise FileNotFoundError(
            "Input directory not found: {}".format(args.in_dir))

    if not os.path.exists(args.qc_report):
        raise FileNotFoundError(
            "QC json file not found: {}".format(args.qc_report))

    if os.path.exists(args.out):
        raise FileExistsError(
            "Output directory already exists: {}".format(args.out))

    with open(args.qc_report, 'r') as qc_report_json:
        qc_report = json.load(qc_report_json)

    files_pass = []
    for bundle_report in qc_report[1]["data"]:
        if bundle_report["status"] == "Pass":
            subject_bundle_fname = str(bundle_report["filename"]).replace(
                ".png",
                ".trk")
            files_pass.append(subject_bundle_fname)

    os.makedirs(args.out)

    for subject_id in os.listdir(args.in_dir):
        os.makedirs(os.path.join(args.out, subject_id))
        bundles_dir = os.path.join(args.in_dir, subject_id, "Clean_Bundles")
        for subject_bundle_fname in os.listdir(bundles_dir):
            possible_fnames = [subject_bundle_fname,
                               subject_bundle_fname.replace("_R_", "_"),
                               subject_bundle_fname.replace("_L_", "_")]

            # Skip file if bundle name is in args.force_remove
            should_skip = False
            for bundle_name_to_remove in args.force_remove:
                for fname in possible_fnames:
                    if bundle_name_to_remove in fname:
                        should_skip = True
            if should_skip:
                continue

            # If file is marked as "Pass", include it
            if any([x in files_pass for x in possible_fnames]):
                bundle_realpath = os.path.realpath(
                    os.path.join(bundles_dir, subject_bundle_fname))
                os.symlink(bundle_realpath, os.path.join(args.out, subject_id,
                                                         subject_bundle_fname))


if __name__ == '__main__':
    main()
