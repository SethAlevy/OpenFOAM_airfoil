import warnings
from pathlib import Path
from typing import List, Literal, Optional

import numpy as np
import pyvista as pv

import utils.utilities as ut
from templates.plot_config import DEFAULT_PLOT_CONFIG, PlotConfig
from utils.logger import SimpleLogger as logger

warnings.filterwarnings("ignore", category=UserWarning, module="pyvista.jupyter")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="pyvista")


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

    if config.pyvista_camera_zoom != 1.0:
        plotter.camera.zoom(config.pyvista_camera_zoom)

    return plotter


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
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    mesh = ut.load_latest_vtm(vtk_dir)

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

    if output_dir is not None:
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

    if output_dir is not None:
        save_pyvista_plot(plotter, output_dir, "velocity_contours", save_formats)
    if show:
        plotter.show()
    plotter.close()


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
    """
    if config is None:
        config = DEFAULT_PLOT_CONFIG

    mesh = ut.load_latest_vtm(vtk_dir)

    p_field = next(
        (name for name in ["p", "p_rgh", "Cp"] if name in mesh.array_names),
        None,
    )

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

    if output_dir is not None:
        save_pyvista_plot(plotter, output_dir, "pressure_contours", save_formats)
    if show:
        plotter.show()
    plotter.close()
