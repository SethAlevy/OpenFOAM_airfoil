import argparse
from pathlib import Path
import pandas as pd
from postprocess.plotting.matplotlib_plots import (
    plot_lift_curve, plot_drag_curve, plot_drag_polar, plot_lift_to_drag_vs_alpha
)
from templates.plot_config import DEFAULT_PLOT_CONFIG
from utils.utilities import read_uiuc_reference_csv


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate polar plots (lift, drag, drag polar, lift-to-drag) from multiple postprocessing summary CSV files."
    )
    parser.add_argument(
        '--summary-csvs',
        nargs='+',
        type=Path,
        required=True,
        help='Paths to postprocessing summary CSV files.'
    )
    parser.add_argument(
        '--reference-csv',
        type=Path,
        required=False,
        help='Optional path to reference CSV file for comparison.'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        required=False,
        default=Path("polar_plots"),
        help='Directory where plots will be saved. Default: ./polar_plots'
    )
    parser.add_argument(
        '--save-formats',
        nargs='+',
        choices=['png', 'pdf'],
        default=['png'],
        help='Output formats for plots. Default: png'
    )
    parser.add_argument(
        '--show',
        action='store_true',
        help='Display plots interactively.'
    )
    return parser.parse_args()


def prepare_series_from_df(df, label, color=None, plot_type="plot"):
    # Use 'Alpha' if present, else 'AngleOfAttack'
    if "Alpha" in df.columns:
        alpha_col = "Alpha"
    elif "AngleOfAttack" in df.columns:
        alpha_col = "AngleOfAttack"
    else:
        raise ValueError("No angle column found (expected 'Alpha' or 'AngleOfAttack').")
    return [
        {"x": df[alpha_col], "y": df["Cl"], "label": f"{label} Cl",
            "color": color, "plot_type": plot_type},
        {"x": df[alpha_col], "y": df["Cd"], "label": f"{label} Cd",
            "color": color, "plot_type": plot_type},
        {"x": df["Cl"], "y": df["Cd"], "label": f"{label} Drag Polar",
            "color": color, "plot_type": plot_type},
        {"x": df[alpha_col], "y": df["Cl"] / df["Cd"],
            "label": f"{label} Cl/Cd", "color": color, "plot_type": plot_type},
    ]


def main():
    args = parse_arguments()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Read and parse all summary CSVs
    all_series = [[], [], [], []]  # For lift, drag, drag polar, lift-to-drag
    colors = ["b", "g", "m", "c", "y", "k", "r"]
    for idx, csv_path in enumerate(args.summary_csvs):
        df = pd.read_csv(csv_path)
        if df is not None:
            label = csv_path.parent.name if csv_path.parent.name else csv_path.stem
            series_list = prepare_series_from_df(
                df, label=label, color=colors[idx % len(colors)])
            for i in range(4):
                all_series[i].append(series_list[i])

    # Step 2: If reference CSV is given, add its data
    if args.reference_csv:
        ref_df = read_uiuc_reference_csv(args.reference_csv)
        ref_series = prepare_series_from_df(
            ref_df, label="Reference", color="orange", plot_type="scatter")
        for i in range(4):
            all_series[i].append(ref_series[i])

    # Step 3: Generate plots
    config = DEFAULT_PLOT_CONFIG
    plot_lift_curve(
        series=all_series[0],
        output_dir=output_dir,
        show=args.show,
        save_formats=args.save_formats,
        config=config
    )
    plot_drag_curve(
        series=all_series[1],
        output_dir=output_dir,
        show=args.show,
        save_formats=args.save_formats,
        config=config
    )
    plot_drag_polar(
        series=all_series[2],
        output_dir=output_dir,
        show=args.show,
        save_formats=args.save_formats,
        config=config
    )
    plot_lift_to_drag_vs_alpha(
        series=all_series[3],
        output_dir=output_dir,
        show=args.show,
        save_formats=args.save_formats,
        config=config
    )
    print(f"Polar plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
