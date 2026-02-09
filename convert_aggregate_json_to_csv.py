#!/usr/bin/env python3

import argparse
import csv
import json
import os

BUNDLES = ['AC',
           'AF_L', 'AF_R',
           'CC_Fr_1', 'CC_Fr_2', 'CC_Oc', 'CC_Pa', 'CC_Pr_Po', 'CC_Te',
           'CG_L', 'CG_R',
           'FAT_L', 'FAT_R',
           'FPT_L', 'FPT_R',
           'FX_L', 'FX_R',
           'ICP_L', 'ICP_R',
           'IFOF_L', 'IFOF_R',
           'ILF_L', 'ILF_R',
           'MCP',
           'MdLF_L', 'MdLF_R',
           'OR_ML_L', 'OR_ML_R',
           'PC',
           'POPT_L', 'POPT_R',
           'SCP_L', 'SCP_R',
           'PYT_L', 'PYT_R',
           'SLF_L', 'SLF_R',
           'UF_L', 'UF_R']
METRICS = ['dice', 'overlap', 'overreach', 'streamlines_count', 'volume',
           'volume_endpoints', 'avg_length', 'std_length', 'min_length',
           'max_length', 'span', 'curl', 'diameter', 'elongation',
           'surface_area', 'end_surface_area_head', 'end_surface_area_tail',
           'irregularity', 'irregularity_of_end_surface_head',
           'irregularity_of_end_surface_tail', 'mean_curvature',
           'fractal_dimension']


def _build_args_parser():
    p = argparse.ArgumentParser(description='',
                                formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('measures_aggregated',
                   help="The final json file 'measures_aggregated.json'.")
    p.add_argument('out_csv',
                   help="Out .csv file")
    p.add_argument('-f', action='store_true',
                   help="Force action. Will overwrite output file if exists.")

    return p


def _read_json(f):
    with open(f, 'r') as result:
        try:
            content = json.load(result)
        except ValueError:
            raise ValueError('File could not be read. json loaded result: {}.'
                             .format(content))
    return content

def _checks(results):
    main_keys = list(results.keys())
    assert main_keys[0] == 'subjects', \
        ('Error in the output file. Expecting a json file similar to:\n'
         '{"subjects": {"subX": {"AF_R": ... }}}')
    subjects = list(results["subjects"])

    return subjects

def _write_to_csv(f, results, subjects):
    """
    Will write an csv file:
    Subject | Bundle | 'dice' | 'overlap' | 'overreach' | ...
    SubjX   | AF_L   |  1.0   |    1.0    |    0.0      | ...
    ...
    """
    with open(f, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Subject", "Bundle"] + METRICS)

        for i, subj in enumerate(subjects):
            nb_vb = 0
            for bundle in BUNDLES:
                if bundle in results[subj]:
                    nb_vb += 1
                    bundle_metrics = []
                    for metric in METRICS:
                        bundle_metrics.append(results[subj][bundle][metric])

                else:
                    bundle_metrics = ['--'] * len(METRICS)
                writer.writerow([subj, bundle] + bundle_metrics)
            print("Subject {}: Found {} valid bundles.".format(subj, nb_vb))

def main():
    parser = _build_args_parser()
    args = parser.parse_args()

    if os.path.isfile(args.out_csv) and not args.f:
        parser.error("Output csv file {} already exists. Can't overwrite."
                     .format(args.out_csv))

    if not os.path.isfile(args.measures_aggregated):
        parser.error("Measures_aggregated file not found: {}"
                     .format(args.measures_aggregated))

    print("Loading results")
    results = _read_json(args.measures_aggregated)
    subjects = _checks(results)
    results = results["subjects"]

    print("Saving to csv information on all {} bundles for {} subjects"
          .format(len(BUNDLES), len(subjects)))
    _write_to_csv(args.out_csv, results, subjects)

    print("Done")

if __name__ == "__main__":
    main()
