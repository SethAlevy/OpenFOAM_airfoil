import pyvista as pv
import matplotlib.pyplot as plt
import utils.utilities as ut
from pathlib import Path
import numpy as np
from templates.plot_config import PlotConfig, DEFAULT_PLOT_CONFIG
from typing import List, Literal, Optional
from utils.logger import SimpleLogger as logger

# Suppress PyVista jupyter backend warnings
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='pyvista.jupyter')
warnings.filterwarnings('ignore', category=DeprecationWarning, module='pyvista')


def _save_matplotlib_plot(
    output_dir: Optional[Path],
    filename: str,
    save_formats: List[Literal['png', 'pdf', 'svg']],
    dpi: int
) -> None:
    """
    Save matplotlib plot in multiple formats.

    Args:
        output_path: Path to save the plot (without extension)
        save_formats: List of formats to save
        dpi: DPI for PNG output
        plot_name: Name of the plot for logging
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
            dpi=dpi if fmt == 'png' else None,
            bbox_inches='tight'
        )
        logger.log(f"Plot saved to: {output_file}")


def _save_pyvista_plot(
    plotter: pv.Plotter,
    output_dir: Optional[Path],
    filename: str,
    save_formats: List[Literal['png', 'html']]
) -> None:
    """Save PyVista plot in multiple formats."""
    if output_dir is None:
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for fmt in save_formats:
        output_file = output_dir / f"{filename}.{fmt}"
        if fmt == 'png':
            plotter.screenshot(str(output_file), transparent_background=False)
            logger.log(f"Plot saved to: {output_file}")
        elif fmt == 'html':
            try:
                plotter.export_html(str(output_file))
                logger.log(f"Plot saved to: {output_file}")
            except ImportError:
                logger.warning(
                    f"Skipping HTML export: trame dependencies not installed. "
                    f"To enable, run: poetry run pip install trame trame-vtk trame-vuetify"
                )
            except Exception as e:
                logger.warning(f"Failed to export HTML: {e}")


def _setup_pyvista_plotter(
    mesh: pv.UnstructuredGrid,
    vtk_dir: Path,
    add_airfoil: bool = True,
    add_mesh: bool = False,
    show: bool = False,
    config: PlotConfig = None
) -> pv.Plotter:
    """Create and configure a PyVista plotter with common settings."""
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    plotter = pv.Plotter(off_screen=not show)

    if add_airfoil:
        airfoil_boundary = ut.load_vtm_boundary(vtk_dir, 'airfoil')
        plotter.add_mesh(
            airfoil_boundary,
            color=config.pyvista_airfoil_color,
            line_width=config.pyvista_airfoil_line_width,
            render_lines_as_tubes=False
        )

    if add_mesh:
        plotter.add_mesh(
            mesh.extract_surface(),
            style='wireframe',
            color=config.pyvista_mesh_color,
            opacity=config.pyvista_mesh_opacity
        )

    # Set XY view and reset camera to show full domain
    plotter.view_xy()
    plotter.reset_camera()  # This ensures the full domain is visible

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
    save_formats: List[Literal['png', 'pdf', 'svg']] = ['png'],
    config: PlotConfig = None,
) -> None:
    """
    Plot the airfoil geometry.

    Args:
        upper_coords: Coordinates of the upper surface.
        lower_coords: Coordinates of the lower surface.
        title: Title of the plot.
        output_dir: Directory to save the plot. Filename will be 'airfoil_geometry'.
        show: Whether to display the plot.
        mean_camber_coords: Coordinates of the mean camber line.
        thickness_coords: Coordinates of the thickness distribution.
        chord: Length of the chord line.
        alpha: Angle of attack in degrees.
        save_formats: List of formats to save. Options: 'png', 'pdf', 'svg'
        config: Plot configuration.
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()

    plt.figure(figsize=config.figsize_airfoil)

    plt.plot(upper_coords[0], upper_coords[1], label='Upper Surface',
             color=config.color_upper_surface, linewidth=config.line_width)
    plt.plot(lower_coords[0], lower_coords[1], label='Lower Surface',
             color=config.color_lower_surface, linewidth=config.line_width)

    if mean_camber_coords is not None:
        plt.plot(mean_camber_coords[0], mean_camber_coords[1], label='Mean Camber Line',
                 color=config.color_camber, linestyle=config.linestyle_camber,
                 linewidth=config.line_width)

    if thickness_coords is not None:
        plt.plot(thickness_coords[0], thickness_coords[1], label='Thickness Distribution',
                 color=config.color_thickness, linestyle=config.linestyle_thickness,
                 linewidth=config.line_width)

    if chord is not None:
        plt.plot([0, np.cos(np.radians(alpha)) * chord],
                 [0, np.sin(np.radians(alpha)) * chord],
                 label='Chord Line', color=config.color_chord,
                 linestyle=config.linestyle_chord, linewidth=config.line_width)

    plt.title(title, fontsize=config.font_size_title,
              fontweight=config.title_fontweight)
    plt.xlabel('Chordwise Location', fontsize=config.font_size_label)
    plt.ylabel('Thickness / Camber', fontsize=config.font_size_label)
    plt.axis('equal')
    plt.grid(True, alpha=config.grid_alpha, linestyle=config.grid_linestyle)
    plt.legend(fontsize=config.font_size_legend)

    if config.tight_layout:
        plt.tight_layout()

    _save_matplotlib_plot(output_dir, 'airfoil_geometry', save_formats, config.dpi)

    if show or output_dir is None:
        plt.show()

    plt.close()


def plot_residuals(
    logs_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal['png', 'pdf', 'svg']] = ['png', 'pdf'],
    config: PlotConfig = None
) -> None:
    """
    Plot residuals from foamLog output files.

    Args:
        logs_dir: Directory containing foamLog output files
        output_dir: Directory to save the plot. Filename will be 'residuals'.
        show: Whether to display the plot
        save_formats: List of formats to save
        config: Plot configuration
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()

    residual_files = {
        'p': 'p_0',
        'Ux': 'Ux_0',
        'Uy': 'Uy_0',
        'k': 'k_0',
        'omega': 'omega_0',
        'epsilon': 'epsilon_0'
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

        ax.semilogy(iterations, residuals, label=label, linewidth=config.line_width,
                    marker='o', markersize=config.marker_size,
                    markevery=max(1, len(iterations) // config.marker_every),
                    color=color)
        plotted_any = True

    if not plotted_any:
        raise FileNotFoundError(f"No residual files found in {logs_dir}")

    ax.set_xlabel('Iteration', fontsize=config.font_size_label)
    ax.set_ylabel('Initial Residual', fontsize=config.font_size_label)
    ax.set_title('Initial Residuals vs Iteration', fontsize=config.font_size_title,
                 fontweight=config.title_fontweight)
    ax.legend(loc='best', fontsize=config.font_size_legend)
    ax.grid(True, which=config.grid_which, linestyle=config.grid_linestyle,
            alpha=config.grid_alpha)
    ax.set_ylim(config.residual_ylim)

    if config.tight_layout:
        plt.tight_layout()

    _save_matplotlib_plot(output_dir, 'residuals', save_formats, config.dpi)

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
    save_formats: List[Literal['png', 'html']] = ['png', 'html'],
    config: PlotConfig = None
) -> None:
    """
    Plot velocity streamlines around airfoil.

    Args:
        vtk_dir: Path to VTK directory
        output_dir: Directory to save the plot. Filename will be 'streamlines'.
        show: Whether to display the plot
        n_streamlines: Number of streamlines to generate
        add_airfoil: Whether to add airfoil boundary
        add_mesh: Whether to show background mesh wireframe
        save_formats: List of formats to save
        config: Plot configuration
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    mesh = ut.load_latest_vtm(vtk_dir)

    if 'U' not in mesh.array_names:
        raise KeyError(f"Velocity field 'U' not found. Available: {mesh.array_names}")

    # Get mesh bounds
    bounds_dict = ut.get_mesh_bounds(mesh)
    x_min, x_max = bounds_dict['x_min'], bounds_dict['x_max']
    y_min, y_max = bounds_dict['y_min'], bounds_dict['y_max']

    # Create seed points
    seed_y = np.linspace(y_min + 0.1 * (y_max - y_min),
                         y_max - 0.1 * (y_max - y_min),
                         n_streamlines)
    seed_x = np.full(n_streamlines, x_min + 0.2 * (x_max - x_min))
    seed_z = np.zeros(n_streamlines)
    seed_points = np.column_stack([seed_x, seed_y, seed_z])

    # Generate streamlines
    streamlines = mesh.streamlines_from_source(
        source=pv.PolyData(seed_points),
        vectors='U',
        terminal_speed=1e-6,
        integration_direction='forward',
        initial_step_length=0.01,
        max_steps=2000,
        compute_vorticity=False
    )

    # Setup plotter - PASS CONFIG HERE
    plotter = _setup_pyvista_plotter(mesh, vtk_dir, add_airfoil, add_mesh, show, config)

    # Add streamlines
    if streamlines.n_points > 0:
        scalar_name = None
        if 'U' in streamlines.array_names:
            streamlines['U_mag'] = ut.get_velocity_magnitude(streamlines)
            scalar_name = 'U_mag'
        elif 'speed' in streamlines.array_names:
            scalar_name = 'speed'

        plotter.add_mesh(
            streamlines,
            scalars=scalar_name,
            cmap=config.pyvista_cmap_velocity,
            line_width=config.pyvista_streamline_line_width,
            render_lines_as_tubes=True,
            scalar_bar_args={
                'title': config.pyvista_scalar_bar_title_velocity} if scalar_name else None
        )
        logger.log(f"Generated {streamlines.n_points} streamline points")

    _save_pyvista_plot(plotter, output_dir, 'streamlines', save_formats)

    if show:
        plotter.show()

    plotter.close()


def plot_velocity_contours(
    vtk_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    add_airfoil: bool = True,
    save_formats: List[Literal['png', 'html']] = ['png', 'html'],
    config: PlotConfig = None
) -> None:
    """
    Plot velocity magnitude contours around airfoil.

    Args:
        vtk_dir: Path to VTK directory
        output_dir: Directory to save the plot. Filename will be 'velocity_contours'.
        show: Display plot
        add_airfoil: Whether to add airfoil boundary
        save_formats: List of formats to save
        config: Plot configuration
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    mesh = ut.load_latest_vtm(vtk_dir)
    mesh['U_mag'] = ut.get_velocity_magnitude(mesh)

    plotter = _setup_pyvista_plotter(mesh, vtk_dir, add_airfoil, False, show, config)

    plotter.add_mesh(
        mesh,
        scalars='U_mag',
        cmap=config.pyvista_cmap_velocity,
        scalar_bar_args={'title': config.pyvista_scalar_bar_title_velocity},
        show_edges=False
    )

    _save_pyvista_plot(plotter, output_dir, 'velocity_contours', save_formats)

    if show:
        plotter.show()

    plotter.close()


def plot_velocity_profiles(
    vtk_dir: Path,
    x_locations: List[float],
    output_dir: Path = None,
    show: bool = False,
    save_formats: List[Literal['png', 'pdf', 'svg']] = ['png', 'pdf'],
    config: PlotConfig = None
) -> None:
    """
    Plot velocity profiles at specified x locations (boundary layer profiles).

    Args:
        vtk_dir: Path to VTK directory
        x_locations: List of x positions to extract profiles
        output_dir: Directory to save the plot. Filename will be 'velocity_profiles'.
        show: Display plot
        save_formats: List of formats to save
        config: Plot configuration
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    config.apply_global_settings()

    mesh = ut.load_latest_vtm(vtk_dir)

    if 'U' not in mesh.array_names:
        raise KeyError(f"Velocity field 'U' not found. Available: {mesh.array_names}")

    fig, ax = plt.subplots(figsize=config.figsize_default)
    colors = plt.cm.viridis(np.linspace(0, 1, len(x_locations)))
    bounds_dict = ut.get_mesh_bounds(mesh)
    y_min, y_max = bounds_dict['y_min'], bounds_dict['y_max']

    for idx, x_loc in enumerate(x_locations):
        # Create sampling line
        line_points = np.zeros((100, 3))
        line_points[:, 0] = x_loc
        line_points[:, 1] = np.linspace(y_min, y_max, 100)

        sampled = pv.PolyData(line_points).sample(mesh)
        U = sampled['U']
        y_coords = sampled.points[:, 1]
        u_mag = np.linalg.norm(U, axis=1) if U.ndim > 1 else U
        valid = ~np.isnan(u_mag)

        ax.plot(u_mag[valid], y_coords[valid], label=f'x = {x_loc:.2f}',
                linewidth=config.line_width, color=colors[idx])

    ax.set_xlabel('Velocity Magnitude [m/s]', fontsize=config.font_size_label)
    ax.set_ylabel('y [m]', fontsize=config.font_size_label)
    ax.set_title('Velocity Profiles at Different x Locations',
                 fontsize=config.font_size_title, fontweight=config.title_fontweight)
    ax.legend(fontsize=config.font_size_legend)
    ax.grid(True, alpha=config.grid_alpha, linestyle=config.grid_linestyle)

    if config.tight_layout:
        plt.tight_layout()

    _save_matplotlib_plot(output_dir, 'velocity_profiles', save_formats, config.dpi)

    if show:
        plt.show()

    plt.close()


def plot_pressure_contours(
    vtk_dir: Path,
    output_dir: Path = None,
    show: bool = False,
    add_airfoil: bool = True,
    save_formats: List[Literal['png', 'html']] = ['png', 'html'],
    config: PlotConfig = None
) -> None:
    """
    Plot pressure contours around airfoil.

    Args:
        vtk_dir: Path to VTK directory
        output_dir: Directory to save the plot. Filename will be 'pressure_contours'.
        show: Display plot
        add_airfoil: Whether to add airfoil boundary
        save_formats: List of formats to save
        config: Plot configuration
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    mesh = ut.load_latest_vtm(vtk_dir)

    # Find pressure field
    p_field = next((name for name in ['p', 'p_rgh', 'Cp']
                   if name in mesh.array_names), None)

    if p_field is None:
        raise KeyError(f"No pressure field found. Available: {mesh.array_names}")

    plotter = _setup_pyvista_plotter(mesh, vtk_dir, add_airfoil, False, show, config)

    plotter.add_mesh(
        mesh,
        scalars=p_field,
        cmap=config.pyvista_cmap_pressure,
        scalar_bar_args={
            'title': f'{p_field} [{config.pyvista_scalar_bar_title_pressure.split("[")[1]}'},
        show_edges=False
    )

    _save_pyvista_plot(plotter, output_dir, 'pressure_contours', save_formats)

    if show:
        plotter.show()

    plotter.close()
