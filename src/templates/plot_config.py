from dataclasses import dataclass
from typing import Dict, Tuple, Optional
import matplotlib.pyplot as plt


@dataclass
class PlotConfig:
    """Global configuration for plot formatting and styling."""

    # Figure settings
    figsize_default: Tuple[int, int] = (12, 8)
    figsize_airfoil: Tuple[int, int] = (10, 5)
    dpi: int = 300

    # Font settings
    font_family: str = 'sans-serif'
    font_size_title: int = 14
    font_size_label: int = 12
    font_size_legend: int = 10
    font_size_tick: int = 10
    title_fontweight: str = 'bold'

    # Line settings
    line_width: int = 2
    marker_size: int = 3
    marker_every: int = 50  # Show marker every N points

    # Grid settings
    grid_alpha: float = 0.6
    grid_linestyle: str = '--'
    grid_which: str = 'both'

    # Color schemes
    colors_default: Tuple[str, ...] = ('b', 'r', 'g', 'm', 'c', 'y', 'k')
    colors_residuals: Dict[str, str] = None

    # Residuals plot settings
    residual_ylim: Tuple[float, float] = (1e-10, 1)

    # Airfoil plot settings
    color_upper_surface: str = 'b'
    color_lower_surface: str = 'r'
    color_camber: str = 'g'
    color_thickness: str = 'm'
    color_chord: str = 'k'
    linestyle_camber: str = '--'
    linestyle_thickness: str = ':'
    linestyle_chord: str = '-'

    # PyVista plot settings
    pyvista_cmap_velocity: str = 'jet'
    pyvista_cmap_pressure: str = 'RdBu_r'
    pyvista_airfoil_color: str = 'black'
    pyvista_airfoil_line_width: int = 4
    pyvista_mesh_color: str = 'lightgray'
    pyvista_mesh_opacity: float = 0.1
    pyvista_camera_zoom: float = 1.0  # Changed from 0.8 to 1.0 (no zoom by default)
    pyvista_streamline_line_width: int = 2
    pyvista_scalar_bar_title_velocity: str = 'Velocity Magnitude [m/s]'
    pyvista_scalar_bar_title_pressure: str = 'Pressure [Pa]'

    # Layout settings
    tight_layout: bool = True

    def __post_init__(self):
        """Initialize default color scheme for residuals if not provided."""
        if self.colors_residuals is None:
            self.colors_residuals = {
                'p': '#1f77b4',      # Blue
                'Ux': '#ff7f0e',     # Orange
                'Uy': '#2ca02c',     # Green
                'k': '#d62728',      # Red
                'omega': '#9467bd',  # Purple
                'epsilon': '#8c564b'  # Brown
            }

    def apply_global_settings(self):
        """Apply global matplotlib settings."""
        plt.rcParams['font.family'] = self.font_family
        plt.rcParams['font.size'] = self.font_size_tick
        plt.rcParams['axes.labelsize'] = self.font_size_label
        plt.rcParams['axes.titlesize'] = self.font_size_title
        plt.rcParams['legend.fontsize'] = self.font_size_legend
        plt.rcParams['xtick.labelsize'] = self.font_size_tick
        plt.rcParams['ytick.labelsize'] = self.font_size_tick
        plt.rcParams['figure.dpi'] = self.dpi
        plt.rcParams['savefig.dpi'] = self.dpi
        plt.rcParams['lines.linewidth'] = self.line_width
        plt.rcParams['lines.markersize'] = self.marker_size
        plt.rcParams['grid.alpha'] = self.grid_alpha
        plt.rcParams['grid.linestyle'] = self.grid_linestyle


# Global default instance
DEFAULT_PLOT_CONFIG = PlotConfig()
