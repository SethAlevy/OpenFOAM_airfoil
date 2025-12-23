import warnings
from pathlib import Path
from typing import List, Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv

import utils.utilities as ut
from templates.plot_config import DEFAULT_PLOT_CONFIG, PlotConfig
from utils.logger import SimpleLogger as logger

# Suppress PyVista jupyter backend warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pyvista.jupyter")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="pyvista")


def save_matplotlib_plot(
    output_dir: Optional[Path],
    filename: str,
    save_formats: List[Literal["png", "pdf"]],
    dpi: int,
) -> None:
    """
    Save the current matplotlib figure in multiple formats.

    Args:
        output_dir (Optional[Path]): Directory where the files will be saved. If None,
            nothing is saved.
        filename (str): Output filename without extension.
        save_formats (List[Literal['png', 'pdf']]): List of file formats to save.
        dpi (int): DPI used for raster formats (PNG). Ignored for vector formats.
    """
    if output_dir is None:
        return

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


def save_pyvista_plot(
    plotter: pv.Plotter,
    output_dir: Optional[Path],
    filename: str,
    save_formats: List[Literal["png", "html"]],
) -> None:
    """
    Save a PyVista plotter view in multiple formats.

    Args:
        plotter (pv.Plotter): Active PyVista plotter.
        output_dir (Optional[Path]): Directory where the files will be saved. If None,
            nothing is saved.
        filename (str): Output filename without extension.
        save_formats (List[Literal['png', 'html']]): List of file formats to save.

    Notes:
        - 'png' saves a static screenshot.
        - 'html' exports an interactive HTML view.
    """
    if output_dir is None:
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for fmt in save_formats:
        output_file = output_dir / f"{filename}.{fmt}"
        if fmt == "png":
            plotter.screenshot(str(output_file), transparent_background=False)
            logger.log(f"Plot saved to: {output_file}")
        elif fmt == "html":
            plotter.export_html(str(output_file))
            logger.log(f"Plot saved to: {output_file}")


def setup_pyvista_plotter(
    mesh: pv.UnstructuredGrid,
    vtk_dir: Path,
    add_airfoil: bool = True,
    add_mesh: bool = False,
    show: bool = False,
    config: PlotConfig = None,
) -> pv.Plotter:
    """
    Create and configure a PyVista plotter with common settings.

    Args:
        mesh (pv.UnstructuredGrid): The mesh used to set camera framing and optional
            mesh overlay.
        vtk_dir (Path): Directory containing VTK/VTU/VTM outputs and boundary data.
        add_airfoil (bool): Whether to overlay the airfoil boundary (patch named
            'airfoil').
        add_mesh (bool): Whether to overlay the surface mesh as a wireframe.
        show (bool): If True, create an on-screen plotter; otherwise use off-screen
            rendering (useful for saving PNGs in headless environments).
        config (PlotConfig): Plot styling configuration. If None, defaults are used.

    Returns:
        pv.Plotter: Configured plotter instance (camera set to XY, full-domain view).
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    plotter = pv.Plotter(off_screen=not show)

    if add_airfoil:
        airfoil_boundary = ut.load_vtm_boundary(vtk_dir, "airfoil")
        plotter.add_mesh(
            airfoil_boundary,
            color=config.pyvista_airfoil_color,
            line_width=config.pyvista_airfoil_line_width,
            render_lines_as_tubes=False,
        )

    if add_mesh:
        plotter.add_mesh(
            mesh.extract_surface(),
            style="wireframe",
            color=config.pyvista_mesh_color,
            opacity=config.pyvista_mesh_opacity,
        )

    plotter.view_xy()
    plotter.reset_camera()

    # Optional: apply zoom after showing full domain
    if config.pyvista_camera_zoom != 1.0:
        plotter.camera.zoom(config.pyvista_camera_zoom)

    return plotter


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

    config.apply_global_settings()

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

    save_matplotlib_plot(output_dir, "airfoil_geometry", save_formats, config.dpi)

    if show or output_dir is None:
        plt.show()

    plt.close()


def plot_residuals(
    logs_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
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

    Raises:
        FileNotFoundError: If no residual files are found in the given directory.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()

    residual_files = {
        "p": "p_0",
        "Ux": "Ux_0",
        "Uy": "Uy_0",
        "k": "k_0",
        "omega": "omega_0",
        "epsilon": "epsilon_0",
    }

    fig, ax = plt.subplots(figsize=config.figsize_default)
    plotted_any = False

    for label, filename in residual_files.items():
        file_path = logs_dir / filename
        if not file_path.exists():
            continue

        data = ut.read_foam_log_file(file_path)
        iterations = data[:, 0] if data.ndim > 1 else np.arange(len(data))
        residuals = data[:, 1] if data.ndim > 1 else data
        color = config.colors_residuals.get(label, config.colors_default[0])

        ax.semilogy(
            iterations,
            residuals,
            label=label,
            linewidth=config.line_width,
            marker="o",
            markersize=config.marker_size,
            markevery=max(1, len(iterations) // config.marker_every),
            color=color,
        )
        plotted_any = True

    if not plotted_any:
        raise FileNotFoundError(f"No residual files found in {logs_dir}")

    ax.set_xlabel("Iteration", fontsize=config.font_size_label)
    ax.set_ylabel("Initial Residual", fontsize=config.font_size_label)
    ax.set_title(
        "Initial Residuals vs Iteration",
        fontsize=config.font_size_title,
        fontweight=config.title_fontweight,
    )
    ax.legend(loc="best", fontsize=config.font_size_legend)
    ax.grid(
        True,
        which=config.grid_which,
        linestyle=config.grid_linestyle,
        alpha=config.grid_alpha,
    )
    ax.set_ylim(config.residual_ylim)

    if config.tight_layout:
        plt.tight_layout()

    save_matplotlib_plot(output_dir, "residuals", save_formats, config.dpi)

    if show:
        plt.show()

    plt.close()


def plot_streamlines(
    vtk_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    n_streamlines: int = 20,
    add_airfoil: bool = True,
    add_mesh: bool = False,
    save_formats: List[Literal["png", "html"]] = ["png", "html"],
    config: PlotConfig = None,
) -> None:
    """
    Plot velocity streamlines around the airfoil.

    Args:
        vtk_dir (Path): Path to the VTK directory created by foamToVTK.
        output_dir (Path): Directory where the plot will be saved. If None, the plot
            is not saved.
        show (bool): Whether to display the plot.
        n_streamlines (int): Number of streamline seed points.
        add_airfoil (bool): Whether to overlay the airfoil boundary.
        add_mesh (bool): Whether to overlay the surface mesh as a wireframe.
        save_formats (List[Literal['png', 'html']]): Output formats to save.
        config (PlotConfig): Plot styling configuration. If None, defaults are used.

    Raises:
        KeyError: If the velocity field 'U' is not found in the VTK dataset.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    mesh = ut.load_latest_vtm(vtk_dir)

    if "U" not in mesh.array_names:
        raise KeyError(
            f"Velocity field 'U' not found. Available: {mesh.array_names}"
        )

    bounds_dict = ut.get_mesh_bounds(mesh)
    x_min, x_max = bounds_dict["x_min"], bounds_dict["x_max"]
    y_min, y_max = bounds_dict["y_min"], bounds_dict["y_max"]

    seed_y = np.linspace(
        y_min + 0.1 * (y_max - y_min),
        y_max - 0.1 * (y_max - y_min),
        n_streamlines,
    )
    seed_x = np.full(n_streamlines, x_min + 0.2 * (x_max - x_min))
    seed_z = np.zeros(n_streamlines)
    seed_points = np.column_stack([seed_x, seed_y, seed_z])

    # Generate streamlines
    streamlines = mesh.streamlines_from_source(
        source=pv.PolyData(seed_points),
        vectors="U",
        terminal_speed=1e-6,
        integration_direction="forward",
        initial_step_length=0.01,
        max_steps=2000,
        compute_vorticity=False,
    )

    plotter = setup_pyvista_plotter(
        mesh=mesh,
        vtk_dir=vtk_dir,
        add_airfoil=add_airfoil,
        add_mesh=add_mesh,
        show=show,
        config=config,
    )

    if streamlines.n_points > 0:
        streamlines["U_mag"] = ut.get_velocity_magnitude(streamlines)
        scalar_name = "U_mag"

        plotter.add_mesh(
            streamlines,
            scalars=scalar_name,
            cmap=config.pyvista_cmap_velocity,
            line_width=config.pyvista_streamline_line_width,
            render_lines_as_tubes=True,
            scalar_bar_args=(
                {"title": config.pyvista_scalar_bar_title_velocity}
                if scalar_name
                else None
            ),
        )
        logger.log(f"Generated {streamlines.n_points} streamline points")

    save_pyvista_plot(plotter, output_dir, "streamlines", save_formats)

    if show:
        plotter.show()

    plotter.close()


def plot_velocity_contours(
    vtk_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    add_airfoil: bool = True,
    save_formats: List[Literal["png", "html"]] = ["png", "html"],
    config: PlotConfig = None,
) -> None:
    """
    Plot velocity magnitude contours around the airfoil.

    Args:
        vtk_dir (Path): Path to the VTK directory created by foamToVTK.
        output_dir (Path): Directory where the plot will be saved. If None, the plot
            is not saved.
        show (bool): Whether to display the plot.
        add_airfoil (bool): Whether to overlay the airfoil boundary.
        save_formats (List[Literal['png', 'html']]): Output formats to save.
        config (PlotConfig): Plot styling configuration. If None, defaults are used.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    mesh = ut.load_latest_vtm(vtk_dir)
    mesh["U_mag"] = ut.get_velocity_magnitude(mesh)

    plotter = setup_pyvista_plotter(
        mesh=mesh,
        vtk_dir=vtk_dir,
        add_airfoil=add_airfoil,
        add_mesh=False,
        show=show,
        config=config,
    )

    plotter.add_mesh(
        mesh,
        scalars="U_mag",
        cmap=config.pyvista_cmap_velocity,
        scalar_bar_args={"title": config.pyvista_scalar_bar_title_velocity},
        show_edges=False,
    )

    save_pyvista_plot(plotter, output_dir, "velocity_contours", save_formats)

    if show:
        plotter.show()

    plotter.close()


def plot_velocity_profiles(
    vtk_dir: Path,
    x_locations: List[float],
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
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

    Raises:
        KeyError: If the velocity field 'U' is not found in the VTK dataset.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()

    mesh = ut.load_latest_vtm(vtk_dir)

    if "U" not in mesh.array_names:
        raise KeyError(
            f"Velocity field 'U' not found. Available: {mesh.array_names}"
        )

    fig, ax = plt.subplots(figsize=config.figsize_default)
    colors = plt.cm.viridis(np.linspace(0, 1, len(x_locations)))
    bounds_dict = ut.get_mesh_bounds(mesh)
    y_min, y_max = bounds_dict["y_min"], bounds_dict["y_max"]

    for idx, x_loc in enumerate(x_locations):
        line_points = np.zeros((100, 3))
        line_points[:, 0] = x_loc
        line_points[:, 1] = np.linspace(y_min, y_max, 100)

        sampled = pv.PolyData(line_points).sample(mesh)
        U = sampled["U"]
        y_coords = sampled.points[:, 1]
        u_mag = np.linalg.norm(U, axis=1) if U.ndim > 1 else U
        valid = ~np.isnan(u_mag)

        ax.plot(
            u_mag[valid],
            y_coords[valid],
            label=f"x = {x_loc:.2f}",
            linewidth=config.line_width,
            color=colors[idx],
        )

    ax.set_xlabel("Velocity Magnitude [m/s]", fontsize=config.font_size_label)
    ax.set_ylabel("y [m]", fontsize=config.font_size_label)
    ax.set_title(
        "Velocity Profiles at Different x Locations",
        fontsize=config.font_size_title,
        fontweight=config.title_fontweight,
    )
    ax.legend(fontsize=config.font_size_legend)
    ax.grid(True, alpha=config.grid_alpha, linestyle=config.grid_linestyle)

    if config.tight_layout:
        plt.tight_layout()

    save_matplotlib_plot(output_dir, "velocity_profiles", save_formats, config.dpi)

    if show:
        plt.show()

    plt.close()


def plot_pressure_contours(
    vtk_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    add_airfoil: bool = True,
    save_formats: List[Literal["png", "html"]] = ["png", "html"],
    config: PlotConfig = None,
) -> None:
    """
    Plot pressure contours around the airfoil.

    Args:
        vtk_dir (Path): Path to the VTK directory created by foamToVTK.
        output_dir (Path): Directory where the plot will be saved. If None, the plot
            is not saved.
        show (bool): Whether to display the plot.
        add_airfoil (bool): Whether to overlay the airfoil boundary.
        save_formats (List[Literal['png', 'html']]): Output formats to save.
        config (PlotConfig): Plot styling configuration. If None, defaults are used.

    Raises:
        KeyError: If no supported pressure field is found in the VTK dataset.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    mesh = ut.load_latest_vtm(vtk_dir)

    p_field = next(
        (name for name in ["p", "p_rgh", "Cp"] if name in mesh.array_names),
        None,
    )

    if p_field is None:
        raise KeyError(f"No pressure field found. Available: {mesh.array_names}")

    plotter = setup_pyvista_plotter(
        mesh=mesh,
        vtk_dir=vtk_dir,
        add_airfoil=add_airfoil,
        add_mesh=False,
        show=show,
        config=config,
    )

    plotter.add_mesh(
        mesh,
        scalars=p_field,
        cmap=config.pyvista_cmap_pressure,
        scalar_bar_args={
            "title": (
                f"{p_field} [{config.pyvista_scalar_bar_title_pressure.split('[')[1]}"
            )
        },
        show_edges=False,
    )

    save_pyvista_plot(plotter, output_dir, "pressure_contours", save_formats)

    if show:
        plotter.show()

    plotter.close()


def plot_force_coefficients(
    case_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
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

    Raises:
        FileNotFoundError: If forceCoeffs output file cannot be found.
        KeyError: If required columns are missing.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()

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

    # Calculate index for last 1/3 of data
    n = len(time)
    zoom_start = n - n // 3

    coeffs = [
        ("Cl", cl, config.colors_default[2] if len(
            config.colors_default) > 2 else None),
        ("Cd", cd, config.colors_default[1] if len(
            config.colors_default) > 1 else None),
        (cm_key, cm, config.colors_default[3] if len(
            config.colors_default) > 3 else None),
    ]

    for name, values, color in coeffs:
        fig, axes = plt.subplots(2, 1, figsize=(
            config.figsize_default[0], config.figsize_default[1] * 1.2), sharex=False)

        # Full plot
        axes[0].plot(
            time,
            values,
            linewidth=config.line_width,
            color=color,
        )
        axes[0].set_ylabel(name, fontsize=config.font_size_label)
        axes[0].set_title(f"{name} (full time history)",
                          fontsize=config.font_size_title)
        axes[0].grid(
            True,
            which=config.grid_which,
            linestyle=config.grid_linestyle,
            alpha=config.grid_alpha,
        )

        # Zoomed plot (last 1/3)
        axes[1].plot(
            time[zoom_start:],
            values[zoom_start:],
            linewidth=config.line_width,
            color=color,
        )
        axes[1].set_ylabel(name, fontsize=config.font_size_label)
        axes[1].set_xlabel("Time [s]", fontsize=config.font_size_label)
        axes[1].set_title(f"{name} (last 1/3, zoomed)", fontsize=config.font_size_title)
        axes[1].grid(
            True,
            which=config.grid_which,
            linestyle=config.grid_linestyle,
            alpha=config.grid_alpha,
        )

        if config.tight_layout:
            plt.tight_layout()

        save_matplotlib_plot(
            output_dir, f"force_coefficients_{name.lower()}", save_formats, config.dpi)

        if show:
            plt.show()

        plt.close()


def plot_cp_distribution_sampled(
    case_dir: Path,
    upper_surface_points: np.ndarray,
    lower_surface_points: np.ndarray,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    airfoil_patch: str = "airfoil",
    config: PlotConfig = None,
) -> None:
    """
    Plot Cp distribution for upper and lower surfaces by sampling Cp from the VTM mesh
    at provided surface coordinates.

    Args:
        case_dir (Path): Path to the OpenFOAM case directory.
        upper_surface_points (np.ndarray): Nx3 array of points along the upper surface.
        lower_surface_points (np.ndarray): Nx3 array of points along the lower surface.
        output_dir (Path): Directory where the plot will be saved. If None, the plot is not saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        airfoil_patch (str): Name of the airfoil patch in the VTK output.
        config (PlotConfig): Plot styling configuration. If None, defaults are used.

    Raises:
        KeyError: If Cp field is not found in the mesh.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()

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

    plt.figure(figsize=config.figsize_default)
    plt.plot(
        upper_s, upper_cp,
        label="Upper Surface",
        color=getattr(config, "color_upper_surface", "b"),
        linewidth=config.line_width,
    )
    plt.plot(
        lower_s, lower_cp,
        label="Lower Surface",
        color=getattr(config, "color_lower_surface", "r"),
        linewidth=config.line_width,
    )
    plt.gca().invert_yaxis()
    plt.xlabel("Surface distance [m]", fontsize=config.font_size_label)
    plt.ylabel("$C_p$", fontsize=config.font_size_label)
    plt.title("Pressure Coefficient Distribution",
              fontsize=config.font_size_title, fontweight=config.title_fontweight)
    plt.legend(fontsize=config.font_size_legend)
    plt.grid(True, alpha=config.grid_alpha, linestyle=config.grid_linestyle)

    if config.tight_layout:
        plt.tight_layout()

    save_matplotlib_plot(output_dir, "cp_distribution", save_formats, config.dpi)

    if show or output_dir is None:
        plt.show()

    plt.close()


def plot_lift_curve(
    csv_path: Path,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
) -> None:
    """
    Plot lift curve (Cl vs Angle of Attack) from a summary CSV of multiple cases.

    Args:
        csv_path (Path): Path to the summary CSV file.
        output_dir (Path): Directory where the plot will be saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()

    df = pd.read_csv(csv_path)
    aoa = df["AngleOfAttack"].astype(float)
    cl = df["Cl"].astype(float)

    plt.figure(figsize=config.figsize_default)
    plt.plot(
        aoa, cl,
        marker="o",
        linestyle="-",
        color=config.colors_default[2] if len(config.colors_default) > 2 else "b",
        linewidth=config.line_width,
    )
    plt.xlabel("Angle of Attack [deg]", fontsize=config.font_size_label)
    plt.ylabel("$C_l$", fontsize=config.font_size_label)
    plt.title("Lift Curve ($C_l$ vs Angle of Attack)",
              fontsize=config.font_size_title, fontweight=config.title_fontweight)
    plt.grid(True, alpha=config.grid_alpha, linestyle=config.grid_linestyle)

    if config.tight_layout:
        plt.tight_layout()

    save_matplotlib_plot(output_dir, "lift_curve", save_formats, config.dpi)

    if show or output_dir is None:
        plt.show()
    plt.close()


def plot_drag_polar(
    csv_path: Path,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal["png", "pdf"]] = ["png", "pdf"],
    config: PlotConfig = None,
) -> None:
    """
    Plot drag polar (Cd vs Cl) from a summary CSV of multiple cases.

    Args:
        csv_path (Path): Path to the summary CSV file.
        output_dir (Path): Directory where the plot will be saved.
        show (bool): Whether to display the plot.
        save_formats (List[Literal['png', 'pdf']]): Output formats to save.
        config (PlotConfig): Plot styling configuration.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()

    df = pd.read_csv(csv_path)
    cl = df["Cl"].astype(float)
    cd = df["Cd"].astype(float)

    plt.figure(figsize=config.figsize_default)
    plt.plot(
        cl, cd,
        marker="o",
        linestyle="-",
        color=config.colors_default[1] if len(config.colors_default) > 1 else "r",
        linewidth=config.line_width,
    )
    plt.xlabel("$C_l$", fontsize=config.font_size_label)
    plt.ylabel("$C_d$", fontsize=config.font_size_label)
    plt.title("Drag Polar ($C_d$ vs $C_l$)", fontsize=config.font_size_title,
              fontweight=config.title_fontweight)
    plt.grid(True, alpha=config.grid_alpha, linestyle=config.grid_linestyle)

    if config.tight_layout:
        plt.tight_layout()

    save_matplotlib_plot(output_dir, "drag_polar", save_formats, config.dpi)

    if show or output_dir is None:
        plt.show()
    plt.close()
