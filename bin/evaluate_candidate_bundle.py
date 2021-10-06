#!/usr/bin/env python3
"""
Evaluate a candidate bundle reconstruction of a single subject, given a
reference ground truth (gold standard) bundle.
"""
import argparse
import json
import os

from dipy.io.streamline import load_tractogram
import numpy as np
from scilpy.tractanalysis.streamlines_metrics import compute_tract_counts_map


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )

    p.add_argument("--in", dest="in_bundle", required=True,
                   help="Input tractogram for the bundle to evaluate. "
                        "Must be .trk")
    p.add_argument("--gs", dest="gs_bundle", required=True,
                   help="Tractogram of the gold standard bundle. "
                        "Must be .trk")
    p.add_argument("--out", required=True,
                   help="Output JSON file to store computed measures.")

    return p


def compute_binary_map(tractogram_fname):
    """Compute the binary map of a tractogram.

    Parameters
    ----------
    tractogram_fname : str
        Path to the tractogram

    Returns
    -------
    binary_map : numpy.ndarray
        A 3D volume, built using the dimensions in the header of the tractogram.

    """
    sft = load_tractogram(tractogram_fname, "same")
    sft.to_vox()
    sft.to_corner()
    streamlines = sft.streamlines
    _, dimensions, _, _ = sft.space_attributes
    if not streamlines:
        pass
    binary_map = compute_tract_counts_map(streamlines, dimensions)
    binary_map[binary_map > 0] = 1
    return binary_map


def compute_voxel_pairwise_measures(bundle_fname, gs_fname):
    """Compute comparison measures between two bundles: Overlap, overreach and
    Dice score. If there is absolutely no overlap, the reconstruction is
    considered a failure; the overlap and dice are set to zero, and the
    overreach is left undefined (None).

    Parameters
    ----------
    bundle_fname : str
        Path to the reconstructed bundle.
    gs_fname : str
        Path to the ground truth (gold standard) bundle.

    Returns
    -------
    measures : dict of (str -> float)
        Dictionary mapping measure names to their value.
    """
    bundle_binary_map = compute_binary_map(bundle_fname)
    gs_binary_map = compute_binary_map(gs_fname)

    bundle_indices = np.where(bundle_binary_map.flatten() > 0)[0]
    gs_indices = np.where(gs_binary_map.flatten() > 0)[0]

    tp = len(np.intersect1d(bundle_indices, gs_indices))
    fp = len(np.setdiff1d(bundle_indices, gs_indices))
    fn = len(np.setdiff1d(gs_indices, bundle_indices))

    if tp == 0:
        overlap = 0.
        overreach = None
        dice = 0.
    else:
        # Overlap / Sensitivity
        overlap = tp / float(tp + fn)
        # Overreach
        overreach = fp / float(tp + fn)
        # Dice
        dice = 2 * tp / float(2 * tp + fp + fn)

    return {"dice": dice, "overlap": overlap, "overreach": overreach}


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    if not os.path.exists(args.in_bundle):
        raise FileNotFoundError(
            "Input tractogram not found: {}".format(args.in_bundle))
    if not os.path.exists(args.gs_bundle):
        raise FileNotFoundError(
            "Gold standard tractogram not found: {}".format(args.gs_bundle))
    if os.path.exists(args.out):
        raise FileExistsError("Output file already exists: {}".format(args.out))

    measures_dict = compute_voxel_pairwise_measures(args.in_bundle,
                                                    args.gs_bundle)
    with open(args.out, 'w') as output_json:
        json.dump(measures_dict, output_json)


if __name__ == '__main__':
    main()
