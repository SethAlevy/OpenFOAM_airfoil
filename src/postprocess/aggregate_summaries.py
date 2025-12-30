import argparse
from pathlib import Path
from utils.logger import SimpleLogger as logger
from utils.utilities import combine_postprocessing_summaries


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Aggregate all postprocessing_summary.csv files from subdirectories into a single combined CSV."
    )
    parser.add_argument(
        '--working-path',
        type=Path,
        required=True,
        help='Path containing case directories with postprocessing_summary.csv files.'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output CSV file path (default: <working_path>/combined_postprocessing_summary.csv)'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.output:
        output_path = args.output
    else:
        output_path = args.working_path / "combined_postprocessing_summary.csv"

    logger.log(f"Aggregating summaries in {args.working_path} ...")
    combine_postprocessing_summaries(args.working_path, output_filename=output_path)
    logger.log(f"Aggregation complete. Output saved to {output_path}")


if __name__ == "__main__":
    main()
