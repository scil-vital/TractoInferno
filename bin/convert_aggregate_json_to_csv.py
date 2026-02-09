#!/usr/bin/env python3

import argparse
import csv
import json
import os

import numpy as np

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
OVERALL_METRICS = ['dice', 'overlap', 'overreach', "total_streamlines"]


def _build_args_parser():
    p = argparse.ArgumentParser(description='',
                                formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('measures_aggregated',
                   help="The final json file 'measures_aggregated.json'.")
    p.add_argument('out_csv_bundles',
                   help="Out .csv file for statistics over individual "
                        "bundles.")
    p.add_argument('out_csv_overall',
                   help="Out .csv file for subject-wise mean statistivs over "
                        "all bundles.")
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
    msg = ('Error in the output file. Expecting a json file similar to:\n'
         '{\n'
         '  "subjects": \n'
         '      {"subX": {"AF_R": ... }},\n'
         '  "scores": {"dice": ...}\n'
         '}')
    assert main_keys[0] == 'subjects', msg
    assert main_keys[1] == 'scores', msg
    subjects = list(results["subjects"].keys())

    return subjects

def _write_to_csv_bundles(f, results, subjects):
    """
    Will write a csv file:
    Subject | Bundle | 'dice' | 'overlap' | 'overreach' | ...
    SubjX   | AF_L   |  1.0   |    1.0    |    0.0      | ...
    ...

    Also computes some overall stats
    """
    dices = []
    overlaps = []
    overreachs = []
    nb_streamlines = []
    overall_stats = {}
    with open(f, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Subject", "Bundle"] + METRICS)

        for subj in subjects:
            nb_found = 0
            nb_ref = 0

            for bundle in BUNDLES:
                if bundle in results[subj]:
                    nb_ref += 1
                    if results[subj][bundle]['streamlines_count'] > 0:
                        nb_found += 1

                    # Overall stats
                    dices.append(results[subj][bundle]['dice'])
                    overlaps.append(results[subj][bundle]['overlap'])
                    overreachs.append(results[subj][bundle]['overreach'])
                    nb_streamlines.append(results[subj][bundle]['streamlines_count'])

                    # Prepare metrics for CSV in order
                    bundle_metrics = []
                    for metric in METRICS:
                        bundle_metrics.append(results[subj][bundle][metric])
                else:
                    bundle_metrics = ['--'] * len(METRICS)
                writer.writerow([subj, bundle] + bundle_metrics)
            print("Subject {}: Found {} valid bundles out of "
                  "{} reference bundle.".format(subj, nb_found, nb_ref))
            overall_stats[subj] = [nb_ref, nb_found, np.mean(dices),
                                   np.mean(overlaps), np.mean(overreachs),
                                   np.sum(nb_streamlines)]

    return overall_stats


def _write_to_csv_overall(f, results, subjects, overall_stats):
    with open(f, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Subject", "Nb reference bundles", "Nb VB"] +
                        OVERALL_METRICS)

        percent_vb = []
        for subj in subjects:
            percent_vb.append(overall_stats[subj][1] / overall_stats[subj][0])
            writer.writerow([subj] + overall_stats[subj])

        writer.writerow(["-----"])
        writer.writerow(["", "", "Percentage of recovered VB"] +
                        OVERALL_METRICS)
        writer.writerow(["TOTAL", "39 possible reference bundles",
                         np.mean(percent_vb) * 100, results['dice'],
                         results['overlap'], results['overreach'],
                         results['total_streamlines']])


def main():
    parser = _build_args_parser()
    args = parser.parse_args()

    if os.path.isfile(args.out_csv_bundles) and not args.f:
        parser.error("Output csv file {} already exists. Can't overwrite."
                     .format(args.out_csv_bundles))

    if os.path.isfile(args.out_csv_overall) and not args.f:
        parser.error("Output csv file {} already exists. Can't overwrite."
                     .format(args.out_csv_overall))

    if not os.path.isfile(args.measures_aggregated):
        parser.error("Measures_aggregated file not found: {}"
                     .format(args.measures_aggregated))

    print("Loading results")
    results = _read_json(args.measures_aggregated)
    subjects = _checks(results)

    print("Saving to csv information on all {} bundles for {} subjects"
          .format(len(BUNDLES), len(subjects)))
    overall_stats = _write_to_csv_bundles(
        args.out_csv_bundles, results["subjects"], subjects)

    print("Saving overall stats")
    _write_to_csv_overall(args.out_csv_overall, results["scores"], subjects,
                          overall_stats)

    print("Done")

if __name__ == "__main__":
    main()
