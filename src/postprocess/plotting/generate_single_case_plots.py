import argparse
from pathlib import Path
from typing import List
from postprocess.plotting.matplotlib_plots import (plot_residuals, plot_velocity_profiles,
                                                   plot_force_coefficients)
from postprocess.plotting.pyvista_plots import (plot_streamlines, plot_velocity_contours,
                                                plot_pressure_contours)
from templates.plot_config import DEFAULT_PLOT_CONFIG
from utils.logger import SimpleLogger


def parse_arguments():
    """Parse command line arguments for plot generation."""
    parser = argparse.ArgumentParser(
        description="Generate visualization plots from OpenFOAM airfoil case results."
    )

    parser.add_argument(
        '--case-dir',
        type=Path,
        required=True,
        help='Path to the OpenFOAM case directory containing VTK and logs folders.'
    )

    parser.add_argument(
        '--output-dir',
        type=Path,
        required=False,
        default=None,
        help='Directory where plots will be saved. Defaults to case-dir/plots.'
    )

    parser.add_argument(
        '--plots',
        nargs='+',
        choices=[
            'residuals', 'streamlines', 'velocity-contours', 'pressure-contours',
            'velocity-profiles', 'force-coeffs', 'all'   # <-- ADD 'force-coeffs'
        ],
        default=['all'],
        help='List of plots to generate. Options: residuals, streamlines, '
             'velocity-contours, pressure-contours, velocity-profiles, force-coeffs, all. '
             'Default: all'
    )

    parser.add_argument(
        '--formats',
        nargs='+',
        choices=['png', 'pdf', 'html'],
        default=['png', 'html'],
        help='Output formats for plots. Options: png, pdf, html. '
             'Note: html only works for PyVista plots (streamlines, contours). '
             'Default: png and html'
    )

    parser.add_argument(
        '--show',
        action='store_true',
        help='Display plots in interactive windows (in addition to saving).'
    )

    parser.add_argument(
        '--n-streamlines',
        type=int,
        default=200,
        help='Number of streamlines to generate for streamline plot. Default: 200'
    )

    parser.add_argument(
        '--velocity-profile-locations',
        nargs='+',
        type=float,
        default=[0.1, 0.5, 1.0],
        help='X-coordinates for velocity profile extraction. Default: 0.1 0.5 1.0'
    )

    parser.add_argument(
        '--add-mesh',
        action='store_true',
        help='Add mesh wireframe to PyVista plots.'
    )

    return parser.parse_args()


def validate_case_directory(case_dir: Path) -> tuple[Path, Path]:
    """
    Validate case directory and return paths to VTK and logs directories.

    Args:
        case_dir: Path to OpenFOAM case directory

    Returns:
        Tuple of (vtk_dir, logs_dir)

    Raises:
        FileNotFoundError: If case directory doesn't exist
    """
    if not case_dir.exists():
        raise FileNotFoundError(f"Case directory does not exist: {case_dir}")

    vtk_dir = case_dir / "VTK"
    logs_dir = case_dir / "logs"

    # Warn about missing directories
    if not vtk_dir.exists():
        SimpleLogger.warning(f"VTK directory not found: {vtk_dir}")

    if not logs_dir.exists():
        SimpleLogger.warning(f"Logs directory not found: {logs_dir}")

    return vtk_dir, logs_dir


def get_save_formats(requested_formats: List[str], plot_type: str) -> List[str]:
    """
    Get appropriate save formats for a given plot type.

    Args:
        requested_formats: Formats requested by user
        plot_type: Type of plot ('matplotlib' or 'pyvista')

    Returns:
        List of valid formats for the plot type
    """
    if plot_type == 'matplotlib':
        valid_formats = ['png', 'pdf']
        return [fmt for fmt in requested_formats if fmt in valid_formats]
    elif plot_type == 'pyvista':
        valid_formats = ['png', 'html']
        return [fmt for fmt in requested_formats if fmt in valid_formats]
    return []


def generate_plots():
    """Main function to generate plots based on command line arguments."""
    args = parse_arguments()

    case_dir = args.case_dir.resolve()
    vtk_dir, logs_dir = validate_case_directory(case_dir)

    output_dir = args.output_dir if args.output_dir else case_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    SimpleLogger.log(f"Generating plots for case: {case_dir}")
    SimpleLogger.log(f"Output directory: {output_dir}")

    config = DEFAULT_PLOT_CONFIG

    # Determine which plots to generate
    plots_to_generate = args.plots
    if 'all' in plots_to_generate:
        plots_to_generate = [
            'residuals', 'streamlines', 'velocity-contours',
            'pressure-contours', 'velocity-profiles', 'force-coeffs'  # <-- ADD HERE
        ]

    if 'residuals' in plots_to_generate:
        if logs_dir.exists():
            SimpleLogger.log("Generating residuals plot...")
            formats = get_save_formats(args.formats, 'matplotlib')
            plot_residuals(
                logs_dir=logs_dir,
                output_dir=output_dir,
                show=args.show,
                save_formats=formats,
                config=config
            )
        else:
            SimpleLogger.warning("Skipping residuals plot (logs directory not found)")

    if not vtk_dir.exists():
        SimpleLogger.warning("Skipping PyVista plots (VTK directory not found)")
        return

    pyvista_formats = get_save_formats(args.formats, 'pyvista')

    if 'streamlines' in plots_to_generate:
        SimpleLogger.log("Generating streamlines plot...")
        plot_streamlines(
            vtk_dir=vtk_dir,
            output_dir=output_dir,
            show=args.show,
            n_streamlines=args.n_streamlines,
            add_airfoil=True,
            add_mesh=args.add_mesh,
            save_formats=pyvista_formats,
            config=config
        )

    if 'velocity-contours' in plots_to_generate:
        SimpleLogger.log("Generating velocity contours plot...")
        plot_velocity_contours(
            vtk_dir=vtk_dir,
            output_dir=output_dir,
            show=args.show,
            add_airfoil=True,
            save_formats=pyvista_formats,
            config=config
        )

    if 'pressure-contours' in plots_to_generate:
        SimpleLogger.log("Generating pressure contours plot...")
        plot_pressure_contours(
            vtk_dir=vtk_dir,
            output_dir=output_dir,
            show=args.show,
            add_airfoil=True,
            save_formats=pyvista_formats,
            config=config
        )
    if 'velocity-profiles' in plots_to_generate:
        SimpleLogger.log("Generating velocity profiles plot...")
        matplotlib_formats = get_save_formats(args.formats, 'matplotlib')
        plot_velocity_profiles(
            vtk_dir=vtk_dir,
            x_locations=args.velocity_profile_locations,
            output_dir=output_dir,
            show=args.show,
            save_formats=matplotlib_formats,
            config=config
        )

    if 'force-coeffs' in plots_to_generate:
        SimpleLogger.log("Generating force coefficients plot...")
        matplotlib_formats = get_save_formats(args.formats, 'matplotlib')
        try:
            plot_force_coefficients(
                case_dir=case_dir,
                output_dir=output_dir,
                show=args.show,
                save_formats=matplotlib_formats,
                config=config
            )
        except Exception as e:
            SimpleLogger.warning(f"Skipping force coefficients plot: {e}")

    SimpleLogger.log(f"\nAll plots saved to: {output_dir}")


if __name__ == "__main__":
    generate_plots()
