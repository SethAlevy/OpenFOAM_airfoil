import numpy as np
import requests
from postprocess.visualizations import plot_airfoil
from templates.airfoil_template import Airfoil
from templates.initial_settings_template import Settings
from scipy.interpolate import interp1d
from pathlib import Path
import utils.utilities as ut


class NACA4(Airfoil):
    def __init__(
            self,
            designation: str = None,
            chord_length: float = None,
            resolution: int = None,
            setup: Settings = None
    ):
        """
        Class to generate NACA 4-digit airfoil. Provides coordinates for the surface
        lines, mean camber, thickness distribution, and angle theta along the chord.
        Is capable to apply angle of attack. The leading edge is always on (0, 0).
        Also allows fast plotting of the airfoil geometry.

        Passed arguments are prioritized, but if not provided, the class will try to
        extract them from the Settings object (initial_settings.json).

        Args:
            designation (str): NACA 4-digit designation, e.g., '2412'
            chord_length (float): Chord length of the airfoil in meters.
            resolution (int): Number of points to discretize the airfoil surface along.
            setup (Settings): Settings object to extract initial parameters from.
        """
        if designation is None:
            self.designation = setup.airfoil_settings.get("Designation")
        else:
            self.designation = designation

        if chord_length is None:
            self.chord_length = setup.airfoil_settings.get("Chord")
        else:
            self.chord_length = chord_length

        if resolution is None:
            self.resolution = setup.airfoil_settings.get("Resolution")
        else:
            self.resolution = resolution
        self._alpha = 0.0

        self.x = np.linspace(0, 1, self.resolution)
        self.x_alpha_zero = self.x.copy()
        self.m, self.p, self.t = self.extract_digits()

    def extract_digits(self) -> tuple[float, float, float]:
        """
        Extract the digits from the NACA designation.

        Returns:
            (m, p, t): maximum camber, location of maximum camber, and maximum
                thickness as fractions of chord length.
        """
        m = int(self.designation[0]) / 100.0  # Maximum camber
        p = int(self.designation[1]) / 10.0    # Location of maximum camber
        t = int(self.designation[2:]) / 100.0   # Maximum thickness

        return m, p, t

    def _mean_camber_line(self, x: np.ndarray, m: float, p: float) -> np.ndarray:
        """
        Calculate the mean camber line.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            m (float): maximum camber
            p (float): location of maximum camber

        Returns:
            np.ndarray: y-coordinates of the mean camber line.
        """
        yc = np.where(
            x < p,
            m / (p ** 2) * (2 * p * x - x ** 2),
            m / ((1 - p) ** 2) * ((1 - 2 * p) + 2 * p * x - x ** 2)
        )
        return yc

    def _thickness_distribution(
            self,
            x: np.ndarray,
            t: float,
            closed: bool = True
    ) -> np.ndarray:
        """
        Calculate the thickness distribution.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            t (float): maximum thickness
            closed (bool): whether the trailing edge is closed or not. Impacts
                the last coefficient in the thickness equation.

        Returns:
            np.ndarray: y-coordinates of the thickness distribution.
        """
        x1 = 0.2969
        x2 = -0.1260
        x3 = -0.3516
        x4 = 0.2843
        x5 = -0.1015 if closed else -0.1036

        yt = 5 * t * (
            x1 * np.sqrt(x)
            + x2 * x
            + x3 * x ** 2
            + x4 * x ** 3
            + x5 * x ** 4
        )
        return yt

    def _mean_camber_derivative(self, x: np.ndarray, m: float, p: float) -> np.ndarray:
        """
        Calculate the derivative of the mean camber line.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            m (float): maximum camber
            p (float): location of maximum camber

        Returns:
            np.ndarray: derivative dyc/dx of the mean camber line.
        """
        dyc_dx = np.where(
            x < p,
            (2 * m / (p ** 2)) * (p - x),
            (2 * m / ((1 - p) ** 2)) * (p - x)
        )
        return dyc_dx

    def _theta(self, dyc_dx: np.ndarray) -> np.ndarray:
        """
        Calculate the angle theta.

        Args:
            dyc_dx (np.ndarray): derivative dyc/dx of the mean camber line.

        Returns:
            np.ndarray: angle theta in radians.
        """
        return np.arctan(dyc_dx)

    def _upper_surface(
            self,
            x: np.ndarray,
            yc: np.ndarray,
            yt: np.ndarray,
            theta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate the upper surface coordinates.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            yc (np.ndarray): y-coordinates of the mean camber line
            yt (np.ndarray): y-coordinates of the thickness distribution
            theta (np.ndarray): angle theta in radians

        Returns:
            tuple[np.ndarray, np.ndarray]: (xu, yu) upper surface coordinates
        """
        xu = x - yt * np.sin(theta)
        yu = yc + yt * np.cos(theta)
        return xu, yu

    def _lower_surface(
            self,
            x: np.ndarray,
            yc: np.ndarray,
            yt: np.ndarray,
            theta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate the lower surface coordinates.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            yc (np.ndarray): y-coordinates of the mean camber line
            yt (np.ndarray): y-coordinates of the thickness distribution
            theta (np.ndarray): angle theta in radians

        Returns:
            tuple[np.ndarray, np.ndarray]: (xl, yl) lower surface coordinates
        """
        xl = x + yt * np.sin(theta)
        yl = yc - yt * np.cos(theta)
        return xl, yl

    def _thickness(self, x: np.ndarray) -> np.ndarray:
        """
        Get the thickness distribution coordinates.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)

        Returns:
            np.ndarray: 2D array of thickness distribution coordinates.
        """
        return self._thickness_distribution(x, self.t)

    def set_angle_of_attack(self, alpha: float) -> None:
        """
        Apply angle of attack to the airfoil. Coordinate system assumes positive angle
        is clockwise. Rotation is done around the leading edge (0, 0).

        Args:
            alpha (float): Angle of attack in degrees.
        """
        self._alpha = alpha

    @property
    def mean_camber_line(self):
        """
        Get the mean camber line coordinates eventually rotated by alpha.
        """
        yc = self._mean_camber_line(self.x_alpha_zero, self.m, self.p)
        return ut.rotate_by_alpha(
            -self.alpha,
            self.x_alpha_zero * self.chord_length,
            yc * self.chord_length
        )

    @property
    def upper_surface(self):
        """
        Get the upper surface coordinates eventually rotated by alpha.
        """
        yc = self._mean_camber_line(self.x_alpha_zero, self.m, self.p)
        yt = self._thickness_distribution(self.x_alpha_zero, self.t)
        dyc_dx = self._mean_camber_derivative(self.x_alpha_zero, self.m, self.p)
        theta = self._theta(dyc_dx)
        xu, yu = self._upper_surface(self.x_alpha_zero, yc, yt, theta)
        return ut.rotate_by_alpha(
            -self.alpha,
            xu * self.chord_length,
            yu * self.chord_length
        )

    @property
    def lower_surface(self):
        """
        Get the lower surface coordinates eventually rotated by alpha.
        """
        yc = self._mean_camber_line(self.x_alpha_zero, self.m, self.p)
        yt = self._thickness_distribution(self.x_alpha_zero, self.t)
        dyc_dx = self._mean_camber_derivative(self.x_alpha_zero, self.m, self.p)
        theta = self._theta(dyc_dx)
        xl, yl = self._lower_surface(self.x_alpha_zero, yc, yt, theta)
        return ut.rotate_by_alpha(
            -self.alpha,
            xl * self.chord_length,
            yl * self.chord_length
        )

    @property
    def theta(self):
        """
        Get the angle theta along the chord.
        """
        dyc_dx = self._mean_camber_derivative(self.x_alpha_zero, self.m, self.p)
        return self._theta(dyc_dx)

    @property
    def thickness(self):
        """
        Get the thickness distribution coordinates eventually rotated by alpha.
        """
        thickness = self._thickness(self.x)
        thickness = np.array([
            self.x * self.chord_length,
            thickness * self.chord_length
        ])
        return ut.rotate_by_alpha(-self.alpha, thickness[0], thickness[1])

    @property
    def chord(self):
        """
        Get the chord length.
        """
        return self.chord_length

    @property
    def alpha(self):
        """
        Get the angle of attack.
        """
        return self._alpha

    def plot(
            self, title: str = "NACA Airfoil",
            save_path: Path = None,
            show: bool = True
    ) -> None:
        """
        Plot the airfoil geometry.

        Args:
            title (str): Title of the plot.
            save_path (Path): Path to save the plot image. If None, the plot is shown.
            show (bool): Whether to display the plot.
        """
        plot_airfoil(
            self.upper_surface,
            self.lower_surface,
            title,
            save_path,
            show,
            self.mean_camber_line,
            self.thickness,
            self.chord,
            -self.alpha,
        )

    def to_stl(self, output_path) -> None:
        """
        Export airfoil geometry to STL format.
        """
        ut.export_airfoil_to_stl_trimesh(
            self.upper_surface,
            self.lower_surface,
            output_path
        )


class UIUCAirfoil:

    def __init__(
            self, designation: str = None,
            chord_length: float = None,
            resolution: int = None,
            setup: Settings = None
    ):
        """
        Class to handle airfoils from the UIUC Airfoil Database. After giving
        the designation, it downloads the .dat file, processes the coordinates,
        resamples and provides surface lines, mean camber line, and thickness
        distribution. Is capable to apply angle of attack. The leading edge is
        always on (0, 0). Also allows fast plotting of the airfoil geometry.

        Args:
            designation (str): Designation of the airfoil, e.g., 'naca2412' or 'mh114'
            chord_length (float): Chord length of the airfoil in meters.
            resolution (int): Number of points to discretize the airfoil surface along.
            setup (Settings): Settings object to extract initial parameters from.
        """
        if designation is None:
            self.designation = setup.airfoil_settings.get("Designation")
        else:
            self.designation = designation

        if chord_length is None:
            self.chord_length = setup.airfoil_settings.get("Chord")
        else:
            self.chord_length = chord_length

        if resolution is None:
            self.resolution = setup.airfoil_settings.get("Resolution")
        else:
            self.resolution = resolution
        self._alpha = 0.0

        self._get_airfoil(self.designation)

    def download_uiuc_airfoil(self, designation: str, save_dir: str = None) -> None:
        """
        Download an airfoil .dat file from the UIUC Airfoil Database.
        Args:
            designation (str): Designation of the airfoil, e.g., 'naca2412' or 'mh114'
            save_dir (str): Directory to save the downloaded file (relative to project
                root). Defaults to input/airfoils.
        """ 
        project_root = Path(__file__).resolve().parents[2]
        if save_dir is None:
            save_path_dir = project_root / "input/airfoils"
        else:
            save_path_dir = project_root / save_dir

        save_path_dir.mkdir(parents=True, exist_ok=True)

        url = f"https://m-selig.ae.illinois.edu/ads/coord/{designation}.dat"
        save_path = save_path_dir / f"{designation}.dat"
        response = requests.get(url)
        if response.status_code == 200:
            with open(save_path, "w") as f:
                f.write(response.text)
            print(f"Downloaded {designation} to {save_path}")
            return save_path
        else:
            print(f"Airfoil '{designation}' not found at UIUC database.")

    def load_airfoil_dat(self, file_path: str) -> np.ndarray:
        """
        Load airfoil coordinates from a .dat file.
        Args:
            file_path (str): Path to the .dat file

        Returns:
            np.ndarray: 2D array containing x and y coordinates
        """

        data = np.loadtxt(file_path, skiprows=1)
        x = data[:, 0]
        y = data[:, 1]
        return np.array([x, y])

    def _extract_upper_lower_lines(
            self,
            coords: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Convert the UIUC dat file to get upper and lower surface coordinates.
        Args:
            coords (np.ndarray): 2xN array with x and y coordinates.
        Returns:
            (upper, lower): Tuple of arrays for upper and lower surfaces.
        """
        x_coords, y_coords = coords

        n_points = int(x_coords[0]) + 1

        x_upper = x_coords[1:n_points]
        y_upper = y_coords[1:n_points]

        x_lower = x_coords[n_points:]
        y_lower = y_coords[n_points:]

        upper = np.array([x_upper, y_upper])
        lower = np.array([x_lower, y_lower])
        return upper, lower

    def _get_airfoil(self, designation: str) -> None:
        """
        Get airfoil coordinates from UIUC database, resample coordinates, calculate
        mean camber line and thickness distribution.
        Args:
            designation (str): Designation of the airfoil.
        """
        file_path = self.download_uiuc_airfoil(designation)
        if not file_path:
            raise ValueError(f"Airfoil '{designation}' not downloaded.")
        coords = self.load_airfoil_dat(file_path)
        upper, lower = self._extract_upper_lower_lines(coords)
        self._upper_line_resampled, self._lower_line_resampled = self._resample_airfoil(
            upper, lower)
        self._mean_camber_line = self._calculate_mean_camber_line()
        self._thickness = self._thickness_distribution()

    def _resample_airfoil(
            self,
            upper_line: np.ndarray,
            lower_line: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Resample airfoil coordinates to have n_points along the chord.
        Args:
            upper_line (np.ndarray): Original array of upper surface coordinates.
            lower_line (np.ndarray): Original array of lower surface coordinates.

        Returns:
            (upper_line_resampled, lower_line_resampled): Arrays of resampled upper and
                lower surface coordinates according to the specified number of points.
        """
        upper_line_resampled = ut.resample_line(
            upper_line[0], upper_line[1], self.resolution)

        lower_line_resampled = ut.resample_line(
            lower_line[0], lower_line[1], self.resolution)

        return upper_line_resampled, lower_line_resampled

    def _calculate_mean_camber_line(self) -> np.ndarray:
        """
        Calculate the mean camber line from upper and lower surfaces lines.

        Returns:
            np.ndarray: Calculated array of mean camber line coordinates.
        """
        x_upper, y_upper = self._upper_line_resampled
        x_lower, y_lower = self._lower_line_resampled

        # Ensure both surfaces have the same x-coordinates for averaging
        x_common = np.linspace(0, 1, self.resolution)
        y_upper_interp = interp1d(x_upper, y_upper, kind='linear',
                                  fill_value="extrapolate")(x_common)
        y_lower_interp = interp1d(x_lower, y_lower, kind='linear',
                                  fill_value="extrapolate")(x_common)

        y_mean = (y_upper_interp + y_lower_interp) / 2.0
        return np.array([x_common, y_mean])

    def _thickness_distribution(self) -> np.ndarray:
        """
        Calculate the thickness distribution from upper and lower surfaces lines.

        Returns:

            np.ndarray: Calculated array of thickness distribution coordinates.
        """
        x_upper, y_upper = self._upper_line_resampled
        x_lower, y_lower = self._lower_line_resampled

        # Ensure both surfaces have the same x-coordinates for thickness calculation
        x_common = np.linspace(0, 1, self.resolution)
        y_upper_interp = interp1d(x_upper, y_upper, kind='linear',
                                  fill_value="extrapolate")(x_common)
        y_lower_interp = interp1d(x_lower, y_lower, kind='linear',
                                  fill_value="extrapolate")(x_common)

        thickness = y_upper_interp - np.abs(y_lower_interp)
        return np.array([x_common, thickness])
    
    def set_angle_of_attack(self, alpha: float) -> None:
        """
        Apply angle of attack to the airfoil. Coordinate system assumes positive angle
        is clockwise. Rotation is done around the leading edge (0, 0).

        Args:
            alpha (float): Angle of attack in degrees.
        """
        self._alpha = alpha

    @property
    def upper_surface(self) -> np.ndarray:
        """Get the upper surface coordinates."""
        upper_line = np.array([
            self._upper_line_resampled[0] * self.chord_length,
            self._upper_line_resampled[1] * self.chord_length
        ])
        return ut.rotate_by_alpha(-self._alpha, upper_line[0], upper_line[1])

    @property
    def lower_surface(self) -> np.ndarray:
        """Get the lower surface coordinates."""
        lower_line = np.array([
            self._lower_line_resampled[0] * self.chord_length,
            self._lower_line_resampled[1] * self.chord_length
        ])
        return ut.rotate_by_alpha(-self._alpha, lower_line[0], lower_line[1])

    @property
    def mean_camber_line(self) -> np.ndarray:
        """Get the mean camber line coordinates."""
        mean_camber_line = np.array([
            self._mean_camber_line[0] * self.chord_length,
            self._mean_camber_line[1] * self.chord_length
        ])
        return ut.rotate_by_alpha(
            -self._alpha,
            mean_camber_line[0],
            mean_camber_line[1]
        )

    @property
    def thickness(self) -> np.ndarray:
        """Get the thickness distribution coordinates."""
        thickness = np.array([
            self._thickness[0] * self.chord_length,
            self._thickness[1] * self.chord_length
        ])
        return ut.rotate_by_alpha(-self._alpha, thickness[0], thickness[1])

    @property
    def chord(self) -> float:
        """Get the chord length."""
        return self.chord_length

    @property
    def alpha(self) -> float:
        """Get the angle of attack."""
        return self._alpha

    def plot(
            self,
            title: str = "UIUC Airfoil",
            save_path: Path = None,
            show: bool = True
    ) -> None:
        """
        Plot the airfoil geometry.

        Args:
            title (str): Title of the plot.
            save_path (Path): Path to save the plot image. If None, the plot is shown.
        """
        plot_airfoil(
            self.upper_surface,
            self.lower_surface,
            title,
            save_path,
            show,
            self.mean_camber_line,
            self.thickness,
            self.chord,
            -self._alpha
        )

    def to_stl(self, output_path) -> None:
        """
        Export airfoil geometry to STL format.
        """
        ut.export_airfoil_to_stl_trimesh(
            self.upper_surface,
            self.lower_surface,
            output_path)
