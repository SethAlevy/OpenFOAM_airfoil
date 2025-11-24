import numpy as np
from scipy.interpolate import interp1d
import trimesh


def rotate_by_alpha(
        alpha: float,
        x: np.ndarray,
        y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Rotate coordinates by angle of attack alpha.

    Args:
        alpha (float): angle of attack in degrees
        x (np.ndarray): x-coordinates
        y (np.ndarray): y-coordinates

    Returns:
        tuple[np.ndarray, np.ndarray]: rotated (x_rotated, y_rotated) coordinates
    """
    alpha_rad = np.radians(alpha)
    x_rot = x * np.cos(alpha_rad) - y * np.sin(alpha_rad)
    y_rot = x * np.sin(alpha_rad) + y * np.cos(alpha_rad)
    return x_rot, y_rot


def resample_line(x, y, n_points):
    """
    Resample a line defined by x and y coordinates to have n_points.

    Args:
        x (np.ndarray): Original x-coordinates.
        y (np.ndarray): Original y-coordinates.
        n_points (int): Number of points to resample to.

    Returns:
        tuple[np.ndarray, np.ndarray]: Resampled x and y coordinates.
    """
    interp_func = interp1d(x, y, kind='linear')
    x_new = np.linspace(np.min(x), np.max(x), n_points)
    y_new = interp_func(x_new)
    return x_new, y_new


def validate_airfoil_designation(designation: str) -> tuple[bool, int]:
    """
    Validate if the airfoil designation is a NACA type and get the number of digits.

    Args:
        designation (str): Airfoil designation string.

    Returns:
        tuple[bool, int]: (is_naca, n_digits). is_naca is True if 'NACA' or 'naca' in name, n_digits is the number of digits after 'NACA' (0 if not NACA).
    """
    if "NACA" in designation or "naca" in designation:
        # Remove 'NACA'/'naca' and count digits
        digits = ''.join(
            filter(str.isdigit, designation.replace("NACA", "").replace("naca", "")))
        return True, len(digits)
    return False, 0


def kinematic_viscosity_air(T: float, rho: float) -> float:
    """
    Calculate kinematic viscosity of air at temperature T [K] and density rho [kg/m^3].
    Uses Sutherland's formula for dynamic viscosity.
    """
    mu0 = 1.716e-5      # Reference dynamic viscosity [kg/(m·s)]
    T0 = 273.15         # Reference temperature [K]
    S = 110.4           # Sutherland constant [K]
    mu = mu0 * ((T / T0) ** 1.5) * (T0 + S) / (T + S)
    nu = mu / rho
    return nu


def international_standard_atmosphere(altitude: float, temperature: float, pressure: float, density: float) -> tuple:
    """
    Calculate temperature [K] and density [kg/m^3] at a given altitude [m] using ISA.
    """
    # Constants
    T0 = 288.15      # Sea level temperature [K]
    p0 = 101325      # Sea level pressure [Pa]
    L = 0.0065       # Temperature lapse rate [K/m]
    g = 9.80665      # Gravity [m/s^2]
    M = 0.0289644    # Molar mass of air [kg/mol]
    R = 8.3144598    # Universal gas constant [J/(mol·K)]
    R_specific = 287.05  # Specific gas constant for dry air [J/(kg·K)]

    # Temperature at altitude
    isa_temperature = T0 - L * altitude

    # Pressure at altitude
    isa_pressure = p0 * (1 - L * altitude / T0) ** (g * M / (R * L))

    # Density at altitude
    isa_density = isa_pressure / (R_specific * isa_temperature)

    if temperature is None:
        temperature = isa_temperature
    if pressure is None:
        pressure = isa_pressure
    if density is None:
        density = isa_density

    return temperature, pressure, density


def export_airfoil_to_stl_trimesh(x, y, output_path, thickness=0.01):
    """
    Export an airfoil to an STL file using trimesh.

    Args:
        x (np.ndarray): x-coordinates of the airfoil
        y (np.ndarray): y-coordinates of the airfoil
        output_path (str): path to the output STL file
        thickness (float): thickness of the extrusion

    Returns:
        None
    """
    # Create a 2D path from airfoil coordinates
    airfoil_2d = np.column_stack((x, y))
    path = trimesh.path.Path2D(entities=[trimesh.path.entities.Line(np.arange(len(x)))],
                               vertices=airfoil_2d)
    # Extrude the 2D path to 3D
    meshes = path.extrude(thickness)
    # If meshes is a list, combine them
    if isinstance(meshes, list):
        combined = trimesh.util.concatenate(meshes)
    else:
        combined = meshes
    combined.export(output_path)
