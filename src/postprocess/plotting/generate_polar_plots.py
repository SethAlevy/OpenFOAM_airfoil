import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from templates.python_template_files.plot_config import DEFAULT_PLOT_CONFIG
from utils.utilities import interpolate_airfoil_coefficients
from postprocess.plotting.matplotlib_plots import (
    plot_lift_curve, plot_drag_curve, plot_drag_polar, plot_lift_to_drag_vs_alpha
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate polar plots (lift, drag, drag polar, lift-to-drag) from \
            multiple postprocessing summary CSV files."
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
        help='Optional path to reference CSV file for interpolating coefficients.'
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


def series_from_df(df, label, color=None, plot_type="plot"):
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


def reference_series_from_interpolation(reference_csv, alpha_values, color="orange"):
    """
    Generate reference series by interpolating coefficients for given alpha values.

    Args:
        reference_csv (Path): Path to reference CSV file.
        alpha_values (np.ndarray): Array of angle of attack values to interpolate.
        color (str): Color for reference plots.

    Returns:
        list: List of 4 series dictionaries (lift, drag, drag polar, lift-to-drag).
    """
    cl_interp = []
    cd_interp = []

    for alpha in alpha_values:
        try:
            coeffs = interpolate_airfoil_coefficients(reference_csv, alpha)
            cl_interp.append(coeffs['Cl'])
            cd_interp.append(coeffs['Cd'])
        except ValueError:
            continue

    cl_interp = np.array(cl_interp)
    cd_interp = np.array(cd_interp)
    alpha_interp = alpha_values[:len(cl_interp)]

    return [
        {"x": alpha_interp, "y": cl_interp, "label": "Reference Cl",
            "color": color, "plot_type": "scatter"},
        {"x": alpha_interp, "y": cd_interp, "label": "Reference Cd",
            "color": color, "plot_type": "scatter"},
        {"x": cl_interp, "y": cd_interp, "label": "Reference Drag Polar",
            "color": color, "plot_type": "scatter"},
        {"x": alpha_interp, "y": cl_interp / cd_interp,
            "label": "Reference Cl/Cd", "color": color, "plot_type": "scatter"},
    ]


def main():
    args = parse_arguments()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    all_series = [[], [], [], []]
    colors = ["b", "g", "m", "c", "y", "k", "r"]
    alpha_values = None

    for idx, csv_path in enumerate(args.summary_csvs):
        df = pd.read_csv(csv_path)
        if df is not None:
            label = csv_path.parent.name or csv_path.stem
            series_list = series_from_df(
                df, label=label, color=colors[idx % len(colors)])
            for i in range(4):
                all_series[i].append(series_list[i])

            if alpha_values is None:
                if "Alpha" in df.columns:
                    alpha_values = df["Alpha"].values
                elif "AngleOfAttack" in df.columns:
                    alpha_values = df["AngleOfAttack"].values

    if args.reference_csv and alpha_values is not None:
        ref_series = reference_series_from_interpolation(
            args.reference_csv, alpha_values, color="orange")
        for i in range(4):
            all_series[i].append(ref_series[i])

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
