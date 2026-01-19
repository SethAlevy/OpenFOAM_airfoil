import numpy as np
from pathlib import Path
from typing import Optional, Tuple
import utils.utilities as ut
import utils.geometry as geo
from utils.logger import SimpleLogger
from templates.python_template_files.airfoil_template import Airfoil
from postprocess.plotting.matplotlib_plots import plot_airfoil


class BaseAirfoil(Airfoil):
    """
    Unified base class for both NACA (analytical) and file-based (e.g., UIUC) airfoils.
    Use source="naca" or source="file" to distinguish.

    Provides coordinates for the surface lines, mean camber, thickness distribution,
    and angle theta along the chord. Is capable to apply angle of attack. The leading
    edge is always on (0, 0). Also allows fast plotting of the airfoil geometry.

    Passed arguments are prioritized, but if not provided, the class will try to
    extract them from the Settings object (initial settings json file).

    Args:
        designation (str): NACA designation or file-based name.
        chord_length (float): Chord length of the airfoil in meters.
        resolution (int): Number of points to discretize the airfoil surface along.
        setup (Settings): Settings object to extract initial parameters from.
    """

    def __init__(
        self,
        source: str = "naca",  # "naca" or "file"
        designation: str = None,
        chord_length: float = None,
        resolution: int = None,
        setup=None
    ):
        airfoil_config = setup.airfoil_settings if setup is not None else {}

        self.source = source
        self.designation = ut.resolve_value(
            designation, airfoil_config, "Designation", designation
        )
        self._chord = ut.resolve_value(
            chord_length, airfoil_config, "Chord", chord_length
        )
        self.resolution = ut.resolve_value(
            resolution, airfoil_config, "Resolution", resolution
        )
        self._alpha = 0.0

        self.x = np.linspace(0, 1.01, self.resolution + 5)
        self.x_alpha_zero = self.x.copy()

        if self.source == "naca":
            self.extract_digits()
        elif self.source == "file":
            self.get_airfoil(self.designation)
        else:
            raise ValueError("Unknown airfoil source type.")

        self.airfoil_details()

    def set_angle_of_attack(self, alpha: float) -> None:
        """
        Apply angle of attack to the airfoil. Coordinate system assumes positive angle
        is clockwise. Rotation is done around the leading edge (0, 0).

        Args:
            alpha (float): Angle of attack in degrees.
        """
        self._alpha = alpha

    @property
    def upper_surface(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the upper surface coordinates eventually rotated by alpha.

        Returns:
            tuple[np.ndarray, np.ndarray]: Two arrays with x and y coordinates.
        """
        if self.source == "naca":
            yc = self._mean_camber_line(self.x_alpha_zero, self.m, self.p)
            yt = self._thickness_distribution(self.x_alpha_zero, self.t)
            dyc_dx = self._mean_camber_derivative(self.x_alpha_zero, self.m, self.p)
            theta = self._theta(dyc_dx)
            xu, yu = self._upper_surface(self.x_alpha_zero, yc, yt, theta)
            xl, yl = self._lower_surface(self.x_alpha_zero, yc, yt, theta)
            mask = np.isclose(self.x_alpha_zero, 1.0, atol=1e-3)
            xu_te, yu_te = xu[mask], yu[mask]
            xl_te, yl_te = xl[mask], yl[mask]
            dists = np.sqrt((xu_te - xl_te)**2 + (yu_te - yl_te)**2)
            idx = np.argmin(dists)
            x_te = (xu_te[idx] + xl_te[idx]) / 2
            y_te = (yu_te[idx] + yl_te[idx]) / 2
            xu[-1], yu[-1] = x_te, y_te
            xl[-1], yl[-1] = x_te, y_te
            xu = xu[:self.resolution]
            yu = yu[:self.resolution]
            return geo.rotate_by_alpha(-self.alpha, xu * self.chord, yu * self.chord)
        elif self.source == "file":
            upper_line = np.array([
                self._upper_line_resampled[0] * self.chord,
                self._upper_line_resampled[1] * self.chord
            ])
            return geo.rotate_by_alpha(-self._alpha, upper_line[0], upper_line[1])

    @property
    def lower_surface(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the lower surface coordinates eventually rotated by alpha.

        Returns:
            tuple[np.ndarray, np.ndarray]: Two arrays with x and y coordinates.
        """
        if self.source == "naca":
            yc = self._mean_camber_line(self.x_alpha_zero, self.m, self.p)
            yt = self._thickness_distribution(self.x_alpha_zero, self.t)
            dyc_dx = self._mean_camber_derivative(self.x_alpha_zero, self.m, self.p)
            theta = self._theta(dyc_dx)
            xu, yu = self._upper_surface(self.x_alpha_zero, yc, yt, theta)
            xl, yl = self._lower_surface(self.x_alpha_zero, yc, yt, theta)
            mask = np.isclose(self.x_alpha_zero, 1.0, atol=1e-3)
            xu_te, yu_te = xu[mask], yu[mask]
            xl_te, yl_te = xl[mask], yl[mask]
            dists = np.sqrt((xu_te - xl_te)**2 + (yu_te - yl_te)**2)
            idx = np.argmin(dists)
            x_te = (xu_te[idx] + xl_te[idx]) / 2
            y_te = (yu_te[idx] + yl_te[idx]) / 2
            xu[-1], yu[-1] = x_te, y_te
            xl[-1], yl[-1] = x_te, y_te
            xl = xl[:self.resolution]
            yl = yl[:self.resolution]
            return geo.rotate_by_alpha(-self.alpha, xl * self.chord, yl * self.chord)
        elif self.source == "file":
            lower_line = np.array([
                self._lower_line_resampled[0] * self.chord,
                self._lower_line_resampled[1] * self.chord
            ])
            return geo.rotate_by_alpha(-self._alpha, lower_line[0], lower_line[1])

    @property
    def mean_camber_line(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the mean camber line coordinates eventually rotated by alpha.

        Returns:
            tuple[np.ndarray, np.ndarray]: Two arrays with x and y coordinates.
        """
        if self.source == "naca":
            yc = self._mean_camber_line(self.x_alpha_zero, self.m, self.p)
            return geo.rotate_by_alpha(
                -self.alpha,
                self.x_alpha_zero * self.chord,
                yc * self.chord
            )
        elif self.source == "file":
            mean_camber_line = np.array([
                self._mean_camber_line[0] * self.chord,
                self._mean_camber_line[1] * self.chord
            ])
            return geo.rotate_by_alpha(
                -self._alpha,
                mean_camber_line[0],
                mean_camber_line[1]
            )

    @property
    def thickness(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the thickness distribution coordinates eventually rotated by alpha.

        Returns:
            tuple[np.ndarray, np.ndarray]: Two arrays with x and y coordinates.
        """
        if self.source == "naca":
            thickness = self._thickness_distribution(self.x, self.t)
            thickness = np.array([
                self.x * self.chord,
                thickness * self.chord
            ])
            return geo.rotate_by_alpha(-self.alpha, thickness[0], thickness[1])
        elif self.source == "file":
            thickness = np.array([
                self._thickness[0] * self.chord,
                self._thickness[1] * self.chord
            ])
            return geo.rotate_by_alpha(-self._alpha, thickness[0], thickness[1])

    @property
    def chord(self) -> float:
        """
        Get the chord length.

        Returns:
            float: chord length.
        """
        return self._chord

    @property
    def alpha(self) -> float:
        """
        Get the angle of attack.

        Returns:
            float: angle of attack in degrees.
        """
        return self._alpha

    def plot(
            self,
            title: str = "Airfoil",
            output_dir: Optional[Path] = None,
            show: bool = True
    ) -> None:
        """
        Plot the airfoil geometry.

        Args:
            title (str): Title of the plot.
            output_dir (Optional[Path]): Directory to save the plot image. If None, the
                plot is shown.
            show (bool): Whether to display the plot.
        """
        plot_airfoil(
            self.upper_surface,
            self.lower_surface,
            title,
            output_dir,
            show,
            self.mean_camber_line,
            self.thickness,
            self.chord,
            -self.alpha,
        )

    def to_stl(
            self,
            output_path: Path,
            thickness: float = 0.002,
            dimension: int = 3
    ) -> None:
        """
        Export airfoil geometry to STL format. It can be 2D or 3D by selecting the
        dimension parameter (in general depends on the mesh requirements).

        Args:
            output_path (Path): Path to save the STL file.
            thickness (float): Thickness for 3D extrusion (only used if dimension=3).
            dimension (int): 2 for 2D STL, 3 for extruded 3D STL.
        """
        xu, yu = self.upper_surface
        xl, yl = self.lower_surface
        if dimension == 2:
            geo.export_airfoil_to_stl_2d(xu, yu, xl, yl, output_path)
        elif dimension == 3:
            geo.export_airfoil_to_stl_3d(xu, yu, xl, yl, output_path, thickness)

    def airfoil_details(self) -> None:
        """
        Simple log with the airfoil details.
        """
        SimpleLogger.log(
            f"Airfoil: {getattr(self, 'designation', getattr(self, 'name', ''))}, "
            f"Chord length: {self.chord}, "
            f"Angle of attack: {self.alpha}"
        )
