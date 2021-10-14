#!/usr/bin/env python
"""
Split a full dataset into sub-datasets (trainset/validset/testset) using
symlinks according to a metadata CSV file that assigns a dataset and a new ID
to every existing ID.
"""
import argparse
import csv
import os


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )

    p.add_argument("--in", dest="in_dir", required=True,
                   help="Input directory with all subjects/bundles.")
    p.add_argument("--metadata", required=True,
                   help="Dataset metadata CSV file. Should contain at least 3 "
                        "columns: original_id,new_id,dataset.")
    p.add_argument("--out", required=True,
                   help="Output directory")

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if not os.path.exists(args.in_dir):
        raise FileNotFoundError(
            "Input directory not found: {}".format(args.in_dir))

    if not os.path.exists(args.metadata):
        raise FileNotFoundError(
            "Metadata file not found: {}".format(args.qc_report))

    if os.path.exists(args.out):
        raise FileExistsError(
            "Output directory already exists: {}".format(args.out))

    metadata_list_dict = []
    with open(args.metadata, 'r') as metadata_csv:
        reader = csv.reader(metadata_csv)
        for i, line in enumerate(reader):
            # Read header
            if i == 0:
                keys = line
            else:
                metadata_list_dict.append({k: v for k, v in zip(keys, line)})

    for subject_dict in metadata_list_dict:
        original_id = subject_dict["original_id"]
        new_id = subject_dict["new_id"]
        dataset = subject_dict["dataset"]

        src_dir_path = os.path.realpath(
            os.path.join(args.in_dir, original_id))
        target_dir_path = os.path.join(args.out, dataset, new_id)

        for subdir in os.listdir(src_dir_path):
            # If directory, symlink all files inside
            if os.path.isdir(os.path.join(src_dir_path, subdir)):
                target_subdir_path = os.path.join(target_dir_path, subdir)
                os.makedirs(target_subdir_path, exist_ok=True)

                for src_fname in os.listdir(os.path.join(src_dir_path, subdir)):
                    src_file_path = os.path.join(src_dir_path, subdir, src_fname)

                    target_fname = src_fname.replace(original_id, new_id)
                    target_file_path = os.path.join(target_subdir_path,
                                                    target_fname)
                    os.symlink(src_file_path, target_file_path)
            # If not a directory, symlink files directly
            else:
                src_fname = subdir
                target_file_path = os.path.join(target_dir_path, src_fname)
                src_file_path = os.path.join(src_dir_path, src_fname)
                os.symlink(src_file_path, target_file_path)


if __name__ == '__main__':
    main()
