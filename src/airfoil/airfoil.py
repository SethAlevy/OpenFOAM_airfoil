import numpy as np
import requests
from scipy.interpolate import interp1d
from pathlib import Path
from utils.logger import SimpleLogger
import utils.utilities as ut

from .airfoil_base import BaseAirfoil


class NACA4(BaseAirfoil):
    """
    Class for NACA 4-digit airfoils. Provides coordinates for the surface lines,
    mean camber, thickness distribution, and angle theta along the chord.
    Is capable to apply angle of attack. The leading edge is always on (0, 0).
    Also allows fast plotting of the airfoil geometry.

    Passed arguments are prioritized, but if not provided, the class will try to
    extract them from the Settings object (initial settings json file).

    Args:
        designation (str): NACA 4-digit designation, e.g., '2412'
        chord_length (float): Chord length of the airfoil in meters.
        resolution (int): Number of points to discretize the airfoil surface along.
        setup (Settings): Settings object to extract initial parameters from.
    """

    def __init__(
        self,
        designation: str = None,
        chord_length: float = None,
        resolution: int = None,
        setup=None
    ):
        super().__init__(
            source="naca",
            designation=designation,
            chord_length=chord_length,
            resolution=resolution,
            setup=setup
        )

    def extract_digits(self):
        """
        Extract the digits from the NACA designation.

        Returns:
            (m, p, t): maximum camber, location of maximum camber, and maximum
                thickness as fractions of chord length.
        """
        m = int(self.designation[0]) / 100.0
        p = int(self.designation[1]) / 10.0
        t = int(self.designation[2:]) / 100.0
        self.m, self.p, self.t = m, p, t

    def _mean_camber_line(self, x, m, p):
        """
        Calculate the mean camber line based on the NACA 4-digit designation.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            m (float): maximum camber
            p (float): location of maximum camber

        Returns:
            np.ndarray: y-coordinates of the mean camber line.
        """
        return np.where(
            x < p,
            m / (p ** 2) * (2 * p * x - x ** 2),
            m / ((1 - p) ** 2) * ((1 - 2 * p) + 2 * p * x - x ** 2)
        )

    def _mean_camber_derivative(self, x, m, p):
        """
        Calculate the derivative of the mean camber line.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            m (float): maximum camber
            p (float): location of maximum camber

        Returns:
            np.ndarray: derivative dyc/dx of the mean camber line.
        """
        return np.where(
            x < p,
            (2 * m / (p ** 2)) * (p - x),
            (2 * m / ((1 - p) ** 2)) * (p - x)
        )

    def _thickness_distribution(self, x, t, closed=True):
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

    def _upper_surface(self, x, yc, yt, theta):
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

    def _lower_surface(self, x, yc, yt, theta):
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

    def _theta(self, dyc_dx):
        """
        Calculate the angle theta which the inverse tangent of the mean camber
        derivative.

        Args:
            dyc_dx (np.ndarray): derivative dyc/dx of the mean camber line.

        Returns:
            np.ndarray: angle theta in radians.
        """
        return np.arctan(dyc_dx)


class NACA5(BaseAirfoil):
    """
    Class for NACA 5-digit airfoils. Provides coordinates for the surface lines,
    mean camber, thickness distribution, and angle theta along the chord.
    Is capable to apply angle of attack. The leading edge is always on (0, 0).
    Also allows fast plotting of the airfoil geometry.

    Passed arguments are prioritized, but if not provided, the class will try to
    extract them from the Settings object (initial settings json file).

    Args:
        designation (str): NACA 5-digit designation, e.g., '23012'
        chord_length (float): Chord length of the airfoil in meters.
        resolution (int): Number of points to discretize the airfoil surface along.
        setup (Settings): Settings object to extract initial parameters from.
    """

    def __init__(
        self,
        designation: str = None,
        chord_length: float = None,
        resolution: int = None,
        setup=None
    ):
        super().__init__(
            source="naca",
            designation=designation,
            chord_length=chord_length,
            resolution=resolution,
            setup=setup
        )

    def extract_digits(self):
        """
        Extract the digits from the NACA 5-digit designation.

        Returns:
            (cl, p, q, t): design lift coefficient, position of max camber,
                camber type, and maximum thickness as fractions of chord length.
        """
        cl = int(self.designation[0]) * 0.15
        p = int(self.designation[1]) / 20.0
        q = int(self.designation[2])
        t = int(self.designation[3:]) / 100.0
        self.cl, self.p, self.q, self.t = cl, p, q, t

    def _mean_camber_line(self, x, cl, p, q):
        """
        Calculate the mean camber line based on the NACA 5-digit designation.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            cl (float): design lift coefficient
            p (float): position of max camber
            q (int): camber type (0 for normal, 1 for reflex)

        Returns:
            np.ndarray: y-coordinates of the mean camber line.
        """
        m = p
        if q == 0:
            k1 = 15.957 * cl
            yc = np.where(
                x < m,
                k1 / 6 * (x**3 - 3 * m * x**2 + m**2 * (3 - m) * x),
                k1 * m**3 / 6 * (1 - x)
            )
        else:
            k2 = 51.99 * cl
            a = 0.2025
            yc = np.where(
                x < m,
                k2 / 6 * (
                    x**3 - 3 * m * x**2 + m**2 * (3 - m) * x + a * (m - x)**3
                ),
                k2 * m**3 / 6 * (1 - x)
            )
        return yc

    def _mean_camber_derivative(self, x, cl, p, q):
        """
        Calculate the derivative of the mean camber line.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            cl (float): design lift coefficient
            p (float): position of max camber
            q (int): camber type

        Returns:
            np.ndarray: derivative dyc/dx of the mean camber line.
        """
        m = p
        if q == 0:
            k1 = 15.957 * cl
            dyc_dx = np.where(
                x < m,
                k1 / 6 * (3 * x**2 - 6 * m * x + m**2 * (3 - m)),
                -k1 * m**3 / 6
            )
        else:
            k2 = 51.99 * cl
            a = 0.2025
            dyc_dx = np.where(
                x < m,
                k2 / 6 * (
                    3 * x**2 - 6 * m * x + m**2 * (3 - m) - 3 * a * (m - x)**2
                ),
                -k2 * m**3 / 6
            )
        return dyc_dx

    def _thickness_distribution(self, x, t, closed=True):
        """
        Calculate the thickness distribution.

        Args:
            x (np.ndarray): x-coordinates along the chord (0 to 1)
            t (float): maximum thickness
            closed (bool): whether the trailing edge is closed or not.

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

    def _upper_surface(self, x, yc, yt, theta):
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

    def _lower_surface(self, x, yc, yt, theta):
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

    def _theta(self, dyc_dx):
        """
        Calculate the angle theta which is the inverse tangent of the mean camber
        derivative.

        Args:
            dyc_dx (np.ndarray): derivative dyc/dx of the mean camber line.

        Returns:
            np.ndarray: angle theta in radians.
        """
        return np.arctan(dyc_dx)


class UIUCAirfoil(BaseAirfoil):
    """
    Class to generate airfoil from UIUC Airfoil Database. Provides coordinates for the
    surface lines, mean camber, thickness distribution, and angle theta along the chord.
    Is capable to apply angle of attack. The leading edge is always on (0, 0).
    Also allows fast plotting of the airfoil geometry.

    Passed arguments are prioritized, but if not provided, the class will try to
    extract them from the Settings object (initial settings json file).

    Args:
        name (str): Name of the airfoil in UIUC database, e.g., 'b737d'
        chord_length (float): Chord length of the airfoil in meters.
        resolution (int): Number of points to discretize the airfoil surface along.
        setup (Settings): Settings object to extract initial parameters from.
    """

    def __init__(
        self,
        name: str = None,
        chord_length: float = None,
        resolution: int = None,
        setup=None
    ):
        super().__init__(
            source="file",
            name=name,
            chord_length=chord_length,
            resolution=resolution,
            setup=setup
        )

    def _get_airfoil(self, designation: str) -> None:
        """
        Get airfoil coordinates from UIUC database, resample coordinates, calculate
        mean camber line and thickness distribution.

        Args:
            str: Designation of the airfoil.
        """
        file_path = self.download_uiuc_airfoil(designation)
        if not file_path:
            raise ValueError(
                f"Airfoil '{designation}' not downloaded."
            )
        coords = self.load_airfoil_dat(file_path)
        upper, lower = self._extract_upper_lower_lines(coords)
        self._upper_line_resampled, self._lower_line_resampled = self._resample_airfoil(
            upper, lower
        )
        self._mean_camber_line = self._calculate_mean_camber_line()
        self._thickness = self._thickness_distribution()

    def download_uiuc_airfoil(self, designation: str, save_dir: str = None) -> Path:
        """
        Download an airfoil .dat file from the UIUC Airfoil Database
        (https://m-selig.ae.illinois.edu/ads/coord_database.html).

        Args:
            designation (str): Designation of the airfoil, e.g., 'naca2412' or 'mh114'
            save_dir (str): Directory to save the downloaded file (relative to project
                root). Defaults to input/airfoils.

        Returns:
            Path: Path to the downloaded file if successful, None otherwise.
        """
        project_root = Path(__file__).resolve().parents[2]
        if save_dir is None:
            save_path_dir = project_root / "input/airfoils"
        else:
            save_path_dir = project_root / save_dir

        save_path_dir.mkdir(parents=True, exist_ok=True)

        url = (
            f"https://m-selig.ae.illinois.edu/ads/coord/{designation}.dat"
        )
        save_path = save_path_dir / f"{designation}.dat"
        response = requests.get(url)
        if response.status_code == 200:
            with open(save_path, "w") as f:
                f.write(response.text)
            SimpleLogger.log(
                f"Downloaded {designation} to {save_path}"
            )
            return save_path
        else:
            SimpleLogger.warning(
                f"Airfoil '{designation}' not found at UIUC database."
            )

    def load_airfoil_dat(self, file_path: Path) -> np.ndarray:
        """
        Load airfoil coordinates from a .dat file.

        Args:
            file_path (Path): Path to the .dat file

        Returns:
            np.ndarray: 2D array containing x and y coordinates
        """
        data = np.loadtxt(file_path, skiprows=1)
        x = data[:, 0]
        y = data[:, 1]
        return np.array([x, y])

    def _extract_upper_lower_lines(self, coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Convert the UIUC dat file to get upper and lower surface coordinates.

        Args:
            coords (np.ndarray): 2D array with x and y coordinates.

        Returns:
            tuple[np.ndarray, np.ndarray]: Tuple of 2D arrays for upper and lower
            surfaces.
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

    def _resample_airfoil(
        self, upper_line: np.ndarray, lower_line: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Resample airfoil coordinates to have n_points along the chord.

        Args:
            upper_line (np.ndarray): Original array of upper surface coordinates.
            lower_line (np.ndarray): Original array of lower surface coordinates.

        Returns:
            tuple[np.ndarray, np.ndarray]: Arrays of resampled upper and lower surface
            coordinates.
        """
        upper_line_resampled = ut.resample_line(
            upper_line[0], upper_line[1], self.resolution
        )
        lower_line_resampled = ut.resample_line(
            lower_line[0], lower_line[1], self.resolution
        )
        return upper_line_resampled, lower_line_resampled

    def _calculate_mean_camber_line(self) -> np.ndarray:
        """
        Calculate the mean camber line from upper and lower surfaces lines.

        Returns:
            np.ndarray: Calculated array of mean camber line coordinates.
        """
        x_upper, y_upper = self._upper_line_resampled
        x_lower, y_lower = self._lower_line_resampled
        x_common = np.linspace(0, 1, self.resolution)
        y_upper_interp = interp1d(
            x_upper, y_upper, kind='linear', fill_value="extrapolate"
        )(x_common)
        y_lower_interp = interp1d(
            x_lower, y_lower, kind='linear', fill_value="extrapolate"
        )(x_common)
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
        x_common = np.linspace(0, 1, self.resolution)
        y_upper_interp = interp1d(
            x_upper, y_upper, kind='linear', fill_value="extrapolate"
        )(x_common)
        y_lower_interp = interp1d(
            x_lower, y_lower, kind='linear', fill_value="extrapolate"
        )(x_common)
        thickness = y_upper_interp - np.abs(y_lower_interp)
        return np.array([x_common, thickness])
