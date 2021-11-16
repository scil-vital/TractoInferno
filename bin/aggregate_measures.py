#!/usr/bin/env python3
"""Aggregate multiple JSON files for subject measures into a single JSON,
while also computing average measures and total streamlines count."""
import argparse
import json
import os

import numpy as np


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )

    p.add_argument("--in", dest="in_measures", required=True, nargs="+",
                   help="Input files for multiple subjects and bundles, "
                        "individual and pairwise. Should be named like: "
                        "{subid}__{bname}_{individual|pairwise}_measures.json")
    p.add_argument("--out", required=True,
                   help="Output JSON file with detailed aggregated measures.")

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    for measure_file in args.in_measures:
        if not os.path.exists(measure_file):
            raise FileNotFoundError(
                "Bundle measure file not found: {}".format(measure_file))
    if os.path.exists(args.out):
        raise FileExistsError(
            "Output file already exists: {}".format(args.out))

    results = {
        "subjects": {},
        "scores": {
            "dice": None,
            "overlap": None,
            "overreach": None,
            "total_streamlines": None
        }
    }
    subject_dict = results["subjects"]
    for measure_file in args.in_measures:
        subid, remaining = measure_file.split("__", 1)

        # Extract bundle name
        bname = remaining.replace("_individual_measures.json", "").replace(
            "_pairwise_measures.json", "")

        if subid not in subject_dict:
            subject_dict[subid] = {}
        sub_results = subject_dict[subid]
        if bname not in subject_dict[subid]:
            # Set default values when no bundle has been reconstruction
            sub_results[bname] = {
                "dice": 0.,
                "overlap": 0.,
                "overreach": None,
                "streamlines_count": 0
            }

        with open(measure_file, 'r') as json_file:
            file_output_dict = json.load(json_file)

        # "Squeeze" single-item lists from output dictionary
        file_output_dict = {k: v[0] if isinstance(v, list) else v for k, v in
                            file_output_dict.items()}
        sub_results[bname].update(file_output_dict)

    # Compute avg scores
    for key in ["dice", "overlap", "overreach"]:
        avg_score = np.mean(
            [bundle_scores[key] for sub_results in results["subjects"].values()
             for bundle_scores in sub_results.values()
             if bundle_scores[key] is not None])
        results["scores"][key] = float(avg_score)

    # Compute total streamlines
    total_streamlines = np.sum(
        [bundle_scores["streamlines_count"] for sub_results in
         results["subjects"].values() for bundle_scores in
         sub_results.values()])
    results["scores"]["total_streamlines"] = int(total_streamlines)

    # Save detailed JSON results
    with open(args.out, 'w') as out_json:
        json.dump(results, out_json)


if __name__ == '__main__':
    main()
