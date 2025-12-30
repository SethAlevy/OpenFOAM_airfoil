import argparse
from pathlib import Path
from utils.utilities import export_postprocessing_summary


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Export postprocessing summary CSV for an OpenFOAM airfoil case."
    )
    parser.add_argument(
        '--case-dir',
        type=Path,
        required=True,
        help='Path to the case directory containing the results.'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output CSV file path (default: <case_dir>/postprocessing_summary.csv)'
    )
    parser.add_argument(
        '--mean-n',
        type=int,
        default=1,
        help='Number of last timesteps to average for coefficients (default: 1, i.e. last value only).'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    export_postprocessing_summary(
        case_dir=args.case_dir,
        output_csv=args.output,
        mean_n=args.mean_n
    )


if __name__ == "__main__":
    main()
