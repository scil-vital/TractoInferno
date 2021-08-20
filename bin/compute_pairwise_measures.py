"""
Evaluate bundle reconstructions of a subject, given ground truth
(gold standard) bundles and candidate bundles.

Compares the candidates against all reference bundles, using filenames to match
bundles.
If a bundle is missing from the candidates, it is considered as a wrong
reconstruction.
If a bundle is missing from the reference, the corresponding candidate is
ignored.
"""
import argparse


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )

    p.add_argument("--in", "in_bundles", required=True,
                   help="Input tractograms for all bundles to evaluate.")
    p.add_argument("--ref", "ref_bundles", required=True,
                   help="Reference tractograms for all bundles to evaluate.")

    return p

def main():
    parser = _build_arg_parser()
    args = parser.parse_args()


if __name__ == '__main__':
    main()
