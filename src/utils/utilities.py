import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d


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


def resample_line(
        x: np.ndarray,
        y: np.ndarray,
        n_points: int
) -> tuple[np.ndarray, np.ndarray]:
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


def split_naca_designation(designation: str) -> tuple[bool, str]:
    """
    Split the airfoil designation is is a NACA type and get digits.

    Args:
        designation (str): Airfoil designation string.

    Returns:
        tuple[bool, str]: (is_naca, digits). is_naca is True if 'NACA' or 'naca' in
        designation, digits is the string of digits after 'NACA' (empty string if not
        NACA).
    """
    if "NACA" in designation or "naca" in designation:
        # Remove 'NACA'/'naca' and count digits
        digits = ''.join(
            filter(str.isdigit, designation.replace("NACA", "").replace("naca", "")))
        return True, digits
    return False, digits


def kinematic_viscosity_air(temperature: float, density: float) -> float:
    """
    Calculate kinematic viscosity of air at temperature T [K] and density rho [kg/m^3].
    Uses Sutherland's formula for dynamic viscosity.

    Args:
        temperature (float): Temperature in Kelvin.
        density (float): Density in kg/m^3.

    Returns:
        float: Kinematic viscosity in m^2/s.
    """
    mu0 = 1.716e-5      # Reference dynamic viscosity [kg/(m·s)]
    T0 = 273.15         # Reference temperature [K]
    S = 110.4           # Sutherland constant [K]
    dynamic_viscosity = mu0 * ((temperature / T0) ** 1.5) * (T0 + S) / (temperature + S)
    kinematic_viscosity = dynamic_viscosity / density
    return kinematic_viscosity


def international_standard_atmosphere(
        altitude: float,
        temperature: float,
        pressure: float,
        density: float
) -> tuple:
    """
    Calculate temperature [K], pressure [Pa], and density [kg/m^3] at a given altitude
    [m] using ISA.

    Args:
        altitude (float): Altitude in meters.
        temperature (float): Temperature in Kelvin. If None, calculated from ISA.
        pressure (float): Pressure in Pascals. If None, calculated from ISA.
        density (float): Density in kg/m^3. If None, calculated from ISA.

    Returns:
        tuple: temperature [K], pressure [Pa], density [kg/m^3]
    """
    # Constants
    T0 = 288.15      # Sea level temperature [K]
    p0 = 101325      # Sea level pressure [Pa]
    L = 0.0065       # Temperature lapse rate [K/m]
    g = 9.80665      # Gravity acceleration [m/s^2]
    M = 0.0289644    # Molar mass of air [kg/mol]
    R = 8.3144598    # Universal gas constant [J/(mol·K)]
    R_specific = 287.05  # Specific gas constant for dry air [J/(kg·K)]

    isa_temperature = T0 - L * altitude
    isa_pressure = p0 * (1 - L * altitude / T0) ** (g * M / (R * L))
    isa_density = isa_pressure / (R_specific * isa_temperature)

    if temperature is None:
        temperature = isa_temperature
    if pressure is None:
        pressure = isa_pressure
    if density is None:
        density = isa_density

    return temperature, pressure, density


def export_airfoil_to_stl_ascii(
        x: np.ndarray,
        y: np.ndarray,
        output_path: Path,
        thickness=0.02
) -> None:
    """
    Export a 2D airfoil contour as a 3D extruded ASCII STL file, centered and extruded
    in both +z and -z directions.

    Args:
        x (np.ndarray): x-coordinates of the closed airfoil contour.
        y (np.ndarray): y-coordinates of the closed airfoil contour.
        output_path (Path): STL file path.
        thickness (float): Total extrusion thickness (extruded equally in +z and -z).
    """
    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)
    z0, z1 = -thickness / 2, thickness / 2

    # Top and bottom vertices
    top = np.column_stack((x, y, np.full_like(x, z1)))
    bottom = np.column_stack((x, y, np.full_like(x, z0)))

    with open(output_path, 'w') as stl:
        stl.write("solid airfoil\n")
        for i in range(n - 1):
            # Side faces (two triangles per quad)
            v0, v1, v2, v3 = top[i], top[i + 1], bottom[i + 1], bottom[i]
            # First triangle
            stl.write("  facet normal 0 0 0\n    outer loop\n")
            stl.write(f"      vertex {v0[0]} {v0[1]} {v0[2]}\n")
            stl.write(f"      vertex {v1[0]} {v1[1]} {v1[2]}\n")
            stl.write(f"      vertex {v2[0]} {v2[1]} {v2[2]}\n")
            stl.write("    endloop\n  endfacet\n")
            # Second triangle
            stl.write("  facet normal 0 0 0\n    outer loop\n")
            stl.write(f"      vertex {v0[0]} {v0[1]} {v0[2]}\n")
            stl.write(f"      vertex {v2[0]} {v2[1]} {v2[2]}\n")
            stl.write(f"      vertex {v3[0]} {v3[1]} {v3[2]}\n")
            stl.write("    endloop\n  endfacet\n")
        stl.write("endsolid airfoil\n")
