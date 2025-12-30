import warnings
from pathlib import Path
from typing import List, Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

import utils.utilities as ut
from templates.plot_config import DEFAULT_PLOT_CONFIG, PlotConfig
from utils.logger import SimpleLogger as logger

warnings.filterwarnings("ignore", category=UserWarning, module="pyvista.jupyter")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="pyvista")


def save_matplotlib_plot(
    output_dir: Path,
    filename: str,
    save_formats: List[Literal["png", "pdf"]],
    dpi: int,
) -> None:
    """
    Save the current matplotlib figure in multiple formats.

    Args:
        output_dir (Path): Directory where the files will be saved.
        filename (str): Output filename without extension.
        save_formats (List[Literal['png', 'pdf']]): List of file formats to save.
        dpi (int): DPI used for raster formats (PNG). Ignored for vector formats.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for fmt in save_formats:
        output_file = output_dir / f"{filename}.{fmt}"
        plt.savefig(
            output_file,
            format=fmt,
            dpi=dpi if fmt == "png" else None,
            bbox_inches="tight",
        )
        logger.log(f"Plot saved to: {output_file}")


def plot_xy_series(
    series: list,
    xlabel: str,
    ylabel: str,
    title: str,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
) -> None:
    """
    Generic XY plot for comparisons.

    Args:
        series: List of dicts, each with keys:
            - "x": x data (array-like)
            - "y": y data (array-like)
            - "label": (optional) label for legend
            - "color": (optional) color
            - "marker": (optional) marker style
            - "plot_type": (optional) "plot" or "scatter" (default "plot")
        xlabel: X axis label.
        ylabel: Y axis label.
        title: Plot title.
        output_dir: Directory to save plot.
        show: Whether to display the plot.
        save_formats: Output formats.
        config: PlotConfig.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()
    plt.figure(figsize=config.figsize_default)

    for s in series:
        x = s["x"]
        y = s["y"]
        label = s.get("label", None)
        color = s.get("color", None)
        marker = s.get("marker", "o")
        s_plot_type = s.get("plot_type", "plot")
        if s_plot_type == "plot":
            plt.plot(x, y, label=label, color=color,
                     marker=marker, linewidth=config.line_width)
        elif s_plot_type == "scatter":
            plt.scatter(x, y, label=label, color=color, marker=marker)

    plt.xlabel(xlabel, fontsize=config.font_size_label)
    plt.ylabel(ylabel, fontsize=config.font_size_label)
    plt.title(title, fontsize=config.font_size_title,
              fontweight=config.title_fontweight)
    plt.grid(True, alpha=config.grid_alpha, linestyle=config.grid_linestyle)
    plt.legend(fontsize=config.font_size_legend)
    if config.tight_layout:
        plt.tight_layout()
    if output_dir is not None:
        save_matplotlib_plot(output_dir, title.replace(
            " ", "_").lower(), save_formats, config.dpi)
    if show or output_dir is None:
        plt.show()
    plt.close()


def plot_airfoil(
    upper_coords: np.ndarray,
    lower_coords: np.ndarray,
    title: str = "Airfoil",
    output_dir: Path = None,
    show: bool = False,
    mean_camber_coords: np.ndarray = None,
    thickness_coords: np.ndarray = None,
    chord: float = None,
    alpha: float = 0,
    save_formats: List[Literal["png", "pdf"]] = ["png"],
    config: PlotConfig = None,
) -> None:
    """
    Plot the airfoil geometry.

    Args:
        upper_coords (np.ndarray): 2xN array of upper surface coordinates.
        lower_coords (np.ndarray): 2xN array of lower surface coordinates.
        title (str): Plot title.
        output_dir (Path): Directory where the plot will be saved. If None, the plot
            is not saved.
        show (bool): Whether to display the plot.
        mean_camber_coords (np.ndarray): Optional 2xN array of mean camber line
            coordinates.
        thickness_coords (np.ndarray): Optional 2xN array of thickness distribution
            coordinates.
        chord (float): Optional chord length for drawing the chord line.
        alpha (float): Angle of attack in degrees (used only for chord-line rotation).
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration. If None, defaults are used.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    plt.figure(figsize=config.figsize_airfoil)
    plt.plot(
        upper_coords[0],
        upper_coords[1],
        label="Upper Surface",
        color=config.color_upper_surface,
        linewidth=config.line_width,
    )
    plt.plot(
        lower_coords[0],
        lower_coords[1],
        label="Lower Surface",
        color=config.color_lower_surface,
        linewidth=config.line_width,
    )
    if mean_camber_coords is not None:
        plt.plot(
            mean_camber_coords[0],
            mean_camber_coords[1],
            label="Mean Camber Line",
            color=config.color_camber,
            linestyle=config.linestyle_camber,
            linewidth=config.line_width,
        )
    if thickness_coords is not None:
        plt.plot(
            thickness_coords[0],
            thickness_coords[1],
            label="Thickness Distribution",
            color=config.color_thickness,
            linestyle=config.linestyle_thickness,
            linewidth=config.line_width,
        )
    if chord is not None:
        plt.plot(
            [0, np.cos(np.radians(alpha)) * chord],
            [0, np.sin(np.radians(alpha)) * chord],
            label="Chord Line",
            color=config.color_chord,
            linestyle=config.linestyle_chord,
            linewidth=config.line_width,
        )
    plt.title(
        title,
        fontsize=config.font_size_title,
        fontweight=config.title_fontweight,
    )
    plt.xlabel("Chordwise Location", fontsize=config.font_size_label)
    plt.ylabel("Thickness / Camber", fontsize=config.font_size_label)
    plt.axis("equal")
    plt.grid(True, alpha=config.grid_alpha, linestyle=config.grid_linestyle)
    plt.legend(fontsize=config.font_size_legend)
    if config.tight_layout:
        plt.tight_layout()
    if show or output_dir is None:
        plt.show()
    if output_dir is not None:
        save_matplotlib_plot(output_dir, "airfoil_geometry", save_formats, config.dpi)
    plt.close()


def plot_residuals(
    logs_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
    plot_type: Literal["plot", "scatter"] = "plot"
) -> None:
    """
    Plot solver residuals from foamLog output files.

    Args:
        logs_dir (Path): Directory containing foamLog output files (e.g. p_0, Ux_0).
        output_dir (Path): Directory where the plot will be saved. If None, the plot
            is not saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration. If None, defaults are used.
        plot_type (Literal["plot", "scatter"]): Plot type for all series.
    """
    residual_files = {
        "p": "p_0",
        "Ux": "Ux_0",
        "Uy": "Uy_0",
        "k": "k_0",
        "omega": "omega_0",
        "epsilon": "epsilon_0",
    }

    series = []
    for label, filename in residual_files.items():
        file_path = logs_dir / filename
        if not file_path.exists():
            continue

        data = ut.read_foam_log_file(file_path)
        iterations = data[:, 0] if data.ndim > 1 else np.arange(len(data))
        residuals = data[:, 1] if data.ndim > 1 else data
        color = config.colors_residuals.get(label, None) if config else None
        series.append({
            "x": iterations,
            "y": residuals,
            "label": label,
            "color": color,
            "marker": "o",
            "plot_type": plot_type
        })

    if not series:
        raise FileNotFoundError(f"No residual files found in {logs_dir}")

    plot_xy_series(
        series=series,
        xlabel="Iteration",
        ylabel="Initial Residual",
        title="Initial Residuals vs Iteration",
        output_dir=output_dir,
        show=show,
        save_formats=save_formats,
        config=config
    )


def plot_velocity_profiles(
    vtk_dir: Path,
    x_locations: List[float],
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
    plot_type: Literal["plot", "scatter"] = "plot"
) -> None:
    """
    Plot velocity profiles at specified x locations (boundary layer profiles).

    Args:
        vtk_dir (Path): Path to the VTK directory created by foamToVTK.
        x_locations (List[float]): List of x positions where velocity profiles will be
            sampled (constant-x lines).
        output_dir (Path): Directory where the plot will be saved. If None, the plot
            is not saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration. If None, defaults are used.
        plot_type (Literal["plot", "scatter"]): Plot type for all series.
    """
    mesh = ut.load_latest_vtm(vtk_dir)

    if "U" not in mesh.array_names:
        raise KeyError(
            f"Velocity field 'U' not found. Available: {mesh.array_names}"
        )

    colors = plt.cm.viridis(np.linspace(0, 1, len(x_locations)))
    bounds_dict = ut.get_mesh_bounds(mesh)
    y_min, y_max = bounds_dict["y_min"], bounds_dict["y_max"]

    series = []
    for idx, x_loc in enumerate(x_locations):
        line_points = np.zeros((100, 3))
        line_points[:, 0] = x_loc
        line_points[:, 1] = np.linspace(y_min, y_max, 100)

        sampled = pv.PolyData(line_points).sample(mesh)
        U = sampled["U"]
        y_coords = sampled.points[:, 1]
        u_mag = np.linalg.norm(U, axis=1) if U.ndim > 1 else U
        valid = ~np.isnan(u_mag)

        series.append({
            "x": u_mag[valid],
            "y": y_coords[valid],
            "label": f"x = {x_loc:.2f}",
            "color": colors[idx],
            "plot_type": plot_type
        })

    plot_xy_series(
        series=series,
        xlabel="Velocity Magnitude [m/s]",
        ylabel="y [m]",
        title="Velocity Profiles at Different x Locations",
        output_dir=output_dir,
        show=show,
        save_formats=save_formats,
        config=config
    )


def plot_force_coefficients(
    case_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
    plot_type: Literal["plot", "scatter"] = "plot"
) -> None:
    """
    Plot force coefficients over simulation time (Cl, Cd, Cm), with zoomed-in views.

    Args:
        case_dir (Path): Path to the OpenFOAM case directory.
        output_dir (Path): Directory where the plot will be saved. If None, the plot
            is not saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration. If None, defaults are used.
        plot_type (Literal["plot", "scatter"]): Plot type for all series.
    """
    file_path = ut.find_latest_force_coeffs_file(case_dir)
    logger.log(f"Reading force coefficients from: {file_path}")

    cols = ut.read_force_coeffs_dat(file_path)

    if "Time" not in cols:
        raise KeyError(f"Missing 'Time' column in: {file_path}")
    if "Cl" not in cols:
        raise KeyError(f"Missing 'Cl' column in: {file_path}")
    if "Cd" not in cols:
        raise KeyError(f"Missing 'Cd' column in: {file_path}")

    time = cols["Time"]
    cl = cols["Cl"]
    cd = cols["Cd"]

    # Prefer CmPitch, then Cm, else first Cm* column.
    cm_key: Optional[str] = None
    if "CmPitch" in cols:
        cm_key = "CmPitch"
    elif "Cm" in cols:
        cm_key = "Cm"
    else:
        for k in cols.keys():
            if k.startswith("Cm"):
                cm_key = k
                break

    if cm_key is None:
        raise KeyError(
            f"No moment coefficient column found in: {file_path}. "
            f"Available columns: {list(cols.keys())}"
        )
    cm = cols[cm_key]

    colors = getattr(config, "colors_default", [None, None, None, None]) if \
        config else [None, None, None, None]
    series = [
        {"x": time, "y": cl, "label": "Cl", "color": colors[2], "plot_type": plot_type},
        {"x": time, "y": cd, "label": "Cd", "color": colors[1], "plot_type": plot_type},
        {"x": time, "y": cm, "label": cm_key, "color": colors[3], "plot_type": plot_type},
    ]

    plot_xy_series(
        series=series,
        xlabel="Time [s]",
        ylabel="Coefficient Value",
        title="Force Coefficients vs Time",
        output_dir=output_dir,
        show=show,
        save_formats=save_formats,
        config=config
    )


def plot_cp_distribution_sampled(
    case_dir: Path,
    upper_surface_points: np.ndarray,
    lower_surface_points: np.ndarray,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
    plot_type: Literal["plot", "scatter"] = "plot"
) -> None:
    """
    Plot Cp distribution for upper and lower surfaces by sampling Cp from the VTM mesh
    at provided surface coordinates.

    Args:
        case_dir (Path): Path to the OpenFOAM case directory.
        upper_surface_points (np.ndarray): Nx3 array of points along the upper surface.
        lower_surface_points (np.ndarray): Nx3 array of points along the lower surface.
        output_dir (Path): Directory where the plot will be saved. If None, the plot is 
        not saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration. If None, defaults are used.
        plot_type (Literal["plot", "scatter"]): Plot type for all series.
    """
    vtk_dir = case_dir / "VTK"
    mesh = ut.load_latest_vtm(vtk_dir)
    if "Cp" not in mesh.array_names:
        raise KeyError("No Cp field found in the VTM mesh.")

    # Sample Cp at upper and lower surface points
    upper_poly = pv.PolyData(upper_surface_points)
    lower_poly = pv.PolyData(lower_surface_points)
    upper_sampled = upper_poly.sample(mesh)
    lower_sampled = lower_poly.sample(mesh)
    upper_cp = upper_sampled["Cp"]
    lower_cp = lower_sampled["Cp"]

    # Use chordwise distance for x-axis (projected onto the surface curve)
    upper_s = np.linalg.norm(upper_surface_points - upper_surface_points[0], axis=1)
    lower_s = np.linalg.norm(lower_surface_points - lower_surface_points[0], axis=1)

    color_upper = getattr(config, "color_upper_surface", "b") if config else "b"
    color_lower = getattr(config, "color_lower_surface", "r") if config else "r"
    series = [
        {"x": upper_s, "y": upper_cp, "label": "Upper Surface", "color": color_upper, "plot_type": plot_type},
        {"x": lower_s, "y": lower_cp, "label": "Lower Surface", "color": color_lower, "plot_type": plot_type},
    ]

    plt.gca().invert_yaxis()
    plot_xy_series(
        series=series,
        xlabel="Surface distance [m]",
        ylabel="$C_p$",
        title="Pressure Coefficient Distribution",
        output_dir=output_dir,
        show=show,
        save_formats=save_formats,
        config=config
    )


def plot_lift_curve(
    series: list,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
) -> None:
    """
    Plot lift curve (Cl vs Angle of Attack) for multiple series.

    Args:
        series: List of dicts, each with keys:
            - "x": Angle of Attack [deg]
            - "y": Cl
            - "label": (optional) label for legend
            - "color": (optional) color
            - "marker": (optional) marker style
            - "plot_type": (optional) "plot" or "scatter"
        output_dir (Path): Directory where the plot will be saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration.
    """
    plot_xy_series(
        series=series,
        xlabel="Angle of Attack [deg]",
        ylabel="$C_l$",
        title="Lift Curve ($C_l$ vs Angle of Attack)",
        output_dir=output_dir,
        show=show,
        save_formats=save_formats,
        config=config
    )


def plot_drag_curve(
    series: list,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
) -> None:
    """
    Plot drag curve (Cd vs Angle of Attack) for multiple series.

    Args:
        series: List of dicts, each with keys:
            - "x": Angle of Attack [deg]
            - "y": Cd
            - "label": (optional) label for legend
            - "color": (optional) color
            - "marker": (optional) marker style
            - "plot_type": (optional) "plot" or "scatter"
        output_dir (Path): Directory where the plot will be saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration.
    """
    plot_xy_series(
        series=series,
        xlabel="Angle of Attack [deg]",
        ylabel="$C_d$",
        title="Drag Curve ($C_d$ vs Angle of Attack)",
        output_dir=output_dir,
        show=show,
        save_formats=save_formats,
        config=config
    )


def plot_drag_polar(
    series: list,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
) -> None:
    """
    Plot drag polar (Cd vs Cl) for multiple series.

    Args:
        series: List of dicts, each with keys:
            - "x": Cl
            - "y": Cd
            - "label": (optional) label for legend
            - "color": (optional) color
            - "marker": (optional) marker style
            - "plot_type": (optional) "plot" or "scatter"
        output_dir (Path): Directory where the plot will be saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration.
    """
    plot_xy_series(
        series=series,
        xlabel="$C_l$",
        ylabel="$C_d$",
        title="Drag Polar ($C_d$ vs $C_l$)",
        output_dir=output_dir,
        show=show,
        save_formats=save_formats,
        config=config
    )


def plot_lift_to_drag_vs_alpha(
    series: list,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
) -> None:
    """
    Plot lift-to-drag ratio (Cl/Cd) vs Angle of Attack for multiple series.

    Args:
        series: List of dicts, each with keys:
            - "x": Angle of Attack [deg]
            - "y": Cl/Cd
            - "label": (optional) label for legend
            - "color": (optional) color
            - "marker": (optional) marker style
            - "plot_type": (optional) "plot" or "scatter"
        output_dir (Path): Directory where the plot will be saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration.
    """
    plot_xy_series(
        series=series,
        xlabel="Angle of Attack [deg]",
        ylabel="$C_l/C_d$",
        title="Lift-to-Drag Ratio ($C_l to C_d$) vs Angle of Attack",
        output_dir=output_dir,
        show=show,
        save_formats=save_formats,
        config=config
    )
