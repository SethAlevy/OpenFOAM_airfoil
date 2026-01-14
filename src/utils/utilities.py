import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
from stl import mesh as np_mesh, Mode
import pyvista as pv
from typing import Optional, List
from utils.logger import SimpleLogger as logger
import pandas as pd
from io import StringIO


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


def resample_by_arc_length(
        x: np.ndarray,
        y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Resample a curve defined by x and y coordinates to have constant arc length
    spacing between points. Returns the same number of points as the input.

    Args:
        x (np.ndarray): Original x-coordinates.
        y (np.ndarray): Original y-coordinates.

    Returns:
        tuple[np.ndarray, np.ndarray]: Resampled x and y coordinates with uniform
        arc length spacing.
    """
    n_points = len(x)

    # Calculate cumulative arc length
    dx = np.diff(x)
    dy = np.diff(y)
    segment_lengths = np.sqrt(dx**2 + dy**2)
    cumulative_length = np.concatenate(([0], np.cumsum(segment_lengths)))

    # Create interpolation functions for x and y as functions of arc length
    interp_x = interp1d(cumulative_length, x, kind='linear')
    interp_y = interp1d(cumulative_length, y, kind='linear')

    # Generate evenly spaced arc length values
    total_length = cumulative_length[-1]
    uniform_arc_lengths = np.linspace(0, total_length, n_points)

    # Interpolate x and y at uniform arc length positions
    x_new = interp_x(uniform_arc_lengths)
    y_new = interp_y(uniform_arc_lengths)

    return x_new, y_new


def export_airfoil_to_stl_3d(
        x_upper: np.ndarray,
        y_upper: np.ndarray,
        x_lower: np.ndarray,
        y_lower: np.ndarray,
        output_path: Path,
        thickness: float = 1e-3,
        resample: bool = True
) -> None:
    """
    Creates and saves a 3D extruded closed airfoil as a numpy-stl ASCII file.
    Optionally resamples coordinates to have constant arc length spacing for better mesh quality.

    Args:
        x_upper (np.ndarray): x-coordinates of the upper surface (LE to TE).
        y_upper (np.ndarray): y-coordinates of the upper surface.
        x_lower (np.ndarray): x-coordinates of the lower surface (TE to LE).
        y_lower (np.ndarray): y-coordinates of the lower surface.
        output_path (Path): STL file path.
        thickness (float): Extrusion thickness in meters (default: 1e-3).
        resample (bool): If True, resample coordinates with uniform arc length spacing (default: True).
    """
    # Resample upper and lower surfaces for uniform arc length spacing
    if resample:
        x_upper, y_upper = resample_by_arc_length(x_upper, y_upper)
        x_lower, y_lower = resample_by_arc_length(x_lower, y_lower)

    x_contour = np.concatenate([x_upper, x_lower[::-1]])
    y_contour = np.concatenate([y_upper, y_lower[::-1]])
    if not (
        np.isclose(x_contour[0], x_contour[-1])
        and np.isclose(y_contour[0], y_contour[-1])
    ):
        x_contour = np.append(x_contour, x_contour[0])
        y_contour = np.append(y_contour, y_contour[0])

    n = len(x_contour)
    z_front = +0.5 * thickness
    z_back = -0.5 * thickness

    front_pts = np.c_[x_contour, y_contour, np.full(n, z_front)]
    back_pts = np.c_[x_contour, y_contour, np.full(n, z_back)]

    faces = []

    faces.extend(
        [front_pts[0], front_pts[i + 1], front_pts[i + 2]]
        for i in range(n - 2)
    )
    faces.extend(
        [back_pts[0], back_pts[i + 2], back_pts[i + 1]] for i in range(n - 2)
    )
    for i in range(n - 1):
        faces.extend(
            (
                [front_pts[i], front_pts[i + 1], back_pts[i + 1]],
                [front_pts[i], back_pts[i + 1], back_pts[i]],
            )
        )
    faces = np.array(faces)
    airfoil_mesh = np_mesh.Mesh(
        np.zeros(faces.shape[0], dtype=np_mesh.Mesh.dtype), name=b'airfoil')
    airfoil_mesh.vectors = faces

    airfoil_mesh.save(str(output_path), mode=Mode.ASCII)


def create_stl_bounding_box(
    x_min: float, x_max: float, y_min: float, y_max: float, z_min: float, z_max: float,
    patch_names: dict, output_dir: Path
) -> None:
    """
    Creates and saves 6 STL files representing the faces of a bounding box.
    Each STL is named according to the patch name.

    Args:
        x_min (float): Minimum x-coordinate of the bounding box.
        x_max (float): Maximum x-coordinate of the bounding box.
        y_min (float): Minimum y-coordinate of the bounding box.
        y_max (float): Maximum y-coordinate of the bounding box.
        z_min (float): Minimum z-coordinate of the bounding box.
        z_max (float): Maximum z-coordinate of the bounding box.
        patch_names (dict): Dictionary mapping patch roles to their names.
        output_dir (Path): Directory to save the STL files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    def make_face(v1, v2, v3, v4, name):
        face_mesh = np_mesh.Mesh(
            np.zeros(2, dtype=np_mesh.Mesh.dtype), name=name.encode())
        face_mesh.vectors[0] = [v1, v2, v3]
        face_mesh.vectors[1] = [v1, v3, v4]
        face_mesh.save(str(output_dir / f"{name}.stl"), mode=Mode.ASCII)

    # Vertices of the box
    p0 = (x_min, y_min, z_min)
    p1 = (x_max, y_min, z_min)
    p2 = (x_max, y_max, z_min)
    p3 = (x_min, y_max, z_min)
    p4 = (x_min, y_min, z_max)
    p5 = (x_max, y_min, z_max)
    p6 = (x_max, y_max, z_max)
    p7 = (x_min, y_max, z_max)

    # Create faces
    make_face(p3, p7, p4, p0, patch_names['inlet'])
    make_face(p1, p5, p6, p2, patch_names['outlet'])
    make_face(p0, p4, p5, p1, patch_names['lowerWall'])
    make_face(p2, p6, p7, p3, patch_names['upperWall'])


def export_airfoil_to_stl_2d(
        x_upper: np.ndarray,
        y_upper: np.ndarray,
        x_lower: np.ndarray,
        y_lower: np.ndarray,
        output_path: Path,
        z_value: float = 0.0
) -> None:
    """
    Export a 2D airfoil contour as a flat ASCII STL file at a single z-value.

    Args:
        x_upper (np.ndarray): x-coordinates of the upper surface.
        y_upper (np.ndarray): y-coordinates of the upper surface.
        x_lower (np.ndarray): x-coordinates of the lower surface.
        y_lower (np.ndarray): y-coordinates of the lower surface.
        output_path (Path): STL file path.
        z_value (float): z-coordinate for all vertices (default: 0.0).
    """
    x_contour = np.concatenate([x_upper, x_lower[::-1]])
    y_contour = np.concatenate([y_upper, y_lower[::-1]])
    n = len(x_contour)

    # Ensure the contour is closed
    if not (np.isclose(x_contour[0], x_contour[-1])
            and np.isclose(y_contour[0], y_contour[-1])):
        x_contour = np.append(x_contour, x_contour[0])
        y_contour = np.append(y_contour, y_contour[0])
        n += 1

    with open(output_path, 'w') as stl:
        stl.write("solid airfoil2d\n")
        # Fan triangulation from the first point
        for i in range(1, n - 1):
            v0 = (x_contour[0], y_contour[0], z_value)
            v1 = (x_contour[i], y_contour[i], z_value)
            v2 = (x_contour[i + 1], y_contour[i + 1], z_value)
            stl.write("  facet normal 0 0 1\n    outer loop\n")
            stl.write(f"      vertex {v0[0]} {v0[1]} {v0[2]}\n")
            stl.write(f"      vertex {v1[0]} {v1[1]} {v1[2]}\n")
            stl.write(f"      vertex {v2[0]} {v2[1]} {v2[2]}\n")
            stl.write("    endloop\n  endfacet\n")
        stl.write("endsolid airfoil2d\n")


def export_domain_to_fms(
    airfoil_x: np.ndarray,
    airfoil_y: np.ndarray,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    z_front: float,
    z_back: float,
    patch_names: dict,
    output_path: Path
) -> None:
    """
    Export the domain around the airfoil as an ASCII FMS file. This type of file
    allows cfMesh to read the patch names for defined faces.

    Args:
        airfoil_x (np.ndarray): x-coordinates of the closed airfoil contour.
        airfoil_y (np.ndarray): y-coordinates of the closed airfoil contour.
        x_min (float): Minimum x-coordinate of the bounding box.
        x_max (float): Maximum x-coordinate of the bounding box.
        y_min (float): Minimum y-coordinate of the bounding box.
        y_max (float): Maximum y-coordinate of the bounding box.
        z_front (float): z-coordinate of the front face.
        z_back (float): z-coordinate of the back face.
        patch_names (dict): Dictionary mapping patch roles to their names.
        output_path (Path): Path to save the FMS file.
    """
    all_patches_with_types = {"airfoil": "wall"}
    all_patches_with_types.update({name: "patch" for name in patch_names.values()})

    num_patches = len(all_patches_with_types)
    patch_to_index = {name: idx for idx,
                      name in enumerate(all_patches_with_types.keys())}

    vertices = []

    num_airfoil_verts = len(airfoil_x)
    vertices.extend(
        (airfoil_x[i], airfoil_y[i], z_front) for i in range(num_airfoil_verts)
    )
    vertices.extend(
        (airfoil_x[i], airfoil_y[i], z_back) for i in range(num_airfoil_verts)
    )

    bbox_front_start = 2 * num_airfoil_verts
    vertices.extend([
        (x_min, y_max, z_front),
        (x_max, y_max, z_front),
        (x_max, y_min, z_front),
        (x_min, y_min, z_front),
    ])

    bbox_back_start = bbox_front_start + 4
    vertices.extend([
        (x_min, y_max, z_back),
        (x_max, y_max, z_back),
        (x_max, y_min, z_back),
        (x_min, y_min, z_back),
    ])

    faces = []

    faces.extend(
        ((0, i, i + 1), patch_to_index["airfoil"])
        for i in range(1, num_airfoil_verts - 1)
    )
    for i in range(1, num_airfoil_verts - 1):
        v0 = num_airfoil_verts  # First vertex on back
        v1 = num_airfoil_verts + i + 1
        v2 = num_airfoil_verts + i
        faces.append(((v0, v1, v2), patch_to_index["airfoil"]))

    for i in range(num_airfoil_verts):
        v0_front = i
        v1_front = (i + 1) % num_airfoil_verts
        v0_back = num_airfoil_verts + i
        v1_back = num_airfoil_verts + (i + 1) % num_airfoil_verts

        faces.append(((v0_front, v1_front, v1_back), patch_to_index["airfoil"]))
        faces.append(((v0_front, v1_back, v0_back), patch_to_index["airfoil"]))

    # Upper wall (top face)
    v0, v1, v2, v3 = bbox_front_start + 0, bbox_front_start + \
        1, bbox_back_start + 1, bbox_back_start + 0
    faces.append(((v0, v1, v2), patch_to_index[patch_names['upperWall']]))
    faces.append(((v0, v2, v3), patch_to_index[patch_names['upperWall']]))

    # Outlet (right face)
    v0, v1, v2, v3 = bbox_front_start + 1, bbox_front_start + \
        2, bbox_back_start + 2, bbox_back_start + 1
    faces.append(((v0, v1, v2), patch_to_index[patch_names['outlet']]))
    faces.append(((v0, v2, v3), patch_to_index[patch_names['outlet']]))

    # Lower wall (bottom face)
    v0, v1, v2, v3 = bbox_front_start + 2, bbox_front_start + \
        3, bbox_back_start + 3, bbox_back_start + 2
    faces.append(((v0, v1, v2), patch_to_index[patch_names['lowerWall']]))
    faces.append(((v0, v2, v3), patch_to_index[patch_names['lowerWall']]))

    # Inlet (left face)
    v0, v1, v2, v3 = bbox_front_start + 3, bbox_front_start + \
        0, bbox_back_start + 0, bbox_back_start + 3
    faces.append(((v0, v1, v2), patch_to_index[patch_names['inlet']]))
    faces.append(((v0, v2, v3), patch_to_index[patch_names['inlet']]))

    with open(output_path, 'w') as f:
        f.write(f"{num_patches}\n(\n\n")
        for name, patch_type in all_patches_with_types.items():
            f.write(f"{name} {patch_type}\n\n")
        f.write(")\n\n\n")

        f.write(f"{len(vertices)}\n(\n")
        for v in vertices:
            f.write(f"({v[0]} {v[1]} {v[2]})\n")
        f.write(")\n\n\n")

        f.write(f"{len(faces)}\n(\n")
        for face_verts, patch_idx in faces:
            f.write(
                f"(({face_verts[0]} {face_verts[1]} {face_verts[2]}) {patch_idx})\n")
        f.write(")\n\n\n")

        edges = set()
        for face_verts, _ in faces:
            for i in range(3):
                v1, v2 = face_verts[i], face_verts[(i + 1) % 3]
                edges.add((min(v1, v2), max(v1, v2)))

        f.write(f"{len(edges)}\n(\n")
        for e in sorted(edges):
            f.write(f"({e[0]} {e[1]})\n")
        f.write(")\n\n")

        f.write("0()\n")
        f.write("0()\n")
        f.write("0()\n")


def read_foam_log_file(file_path: Path) -> np.ndarray:
    """
    Read a foamLog output file and return the data as numpy array.

    Args:
        file_path: Path to the log file

    Returns:
        numpy array with the data

    Raises:
        FileNotFoundError: If the log file does not exist
        ValueError: If the file cannot be parsed or is empty
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Log file not found: {file_path}")

    data = np.loadtxt(file_path, comments='#')
    return data


def load_latest_vtm(vtk_dir: Path) -> pv.UnstructuredGrid:
    """
    Load the latest VTM file from OpenFOAM VTK output directory.

    Args:
        vtk_dir: Path to VTK directory (e.g., case/VTK/)

    Returns:
        PyVista mesh object

    Raises:
        FileNotFoundError: If no VTM files found
        ValueError: If VTM file contains no valid blocks
    """
    vtm_files = sorted(vtk_dir.glob("*_*.vtm"))

    if not vtm_files:
        raise FileNotFoundError(f"No VTM files found in {vtk_dir}")

    latest_file = vtm_files[-1]
    logger.log(f"Loading VTM file: {latest_file.name}")

    multiblock = pv.read(latest_file)
    logger.log(f"Multi-block dataset loaded with {multiblock.n_blocks} blocks")

    # Extract internal mesh
    if 'internal' in multiblock.keys():
        mesh = multiblock['internal']
        logger.log(f"Internal mesh: {mesh.n_cells} cells, {mesh.n_points} points")
    else:
        # Use first non-empty block
        mesh = next((multiblock[i] for i in range(multiblock.n_blocks)
                     if multiblock[i] is not None and multiblock[i].n_cells > 0), None)

        if mesh is None:
            raise ValueError(f"No valid mesh blocks found in {latest_file}")

        logger.log(f"Using mesh with {mesh.n_cells} cells")

    logger.log(f"Available fields: {mesh.array_names}")
    return mesh


def load_vtm_boundary(vtk_dir: Path, patch_name: str = 'airfoil') -> pv.PolyData:
    """
    Load a specific boundary patch from OpenFOAM VTM output.

    Args:
        vtk_dir: Path to VTK directory
        patch_name: Name of boundary patch (e.g., 'airfoil', 'inlet')

    Returns:
        PyVista PolyData of the boundary patch

    Raises:
        FileNotFoundError: If VTM files or boundary file not found
    """
    vtm_files = sorted(vtk_dir.glob("*_*.vtm"))

    if not vtm_files:
        raise FileNotFoundError(f"No VTM files found in {vtk_dir}")

    latest_file = vtm_files[-1]
    time_folder = latest_file.stem
    boundary_file = vtk_dir / time_folder / 'boundary' / f'{patch_name}.vtp'

    if not boundary_file.exists():
        boundary_dir = vtk_dir / time_folder / 'boundary'
        available = [f.stem for f in boundary_dir.glob(
            '*.vtp')] if boundary_dir.exists() else []
        raise FileNotFoundError(
            f"Boundary '{patch_name}' not found. Available: {available}"
        )

    boundary = pv.read(boundary_file)
    logger.log(
        f"Loaded boundary '{patch_name}': {boundary.n_points}"
        f"points, {boundary.n_cells} cells")
    return boundary


def get_velocity_magnitude(mesh: pv.UnstructuredGrid, vel_field: str = 'U') -> np.ndarray:
    """
    Calculate velocity magnitude from velocity vector field.

    Args:
        mesh: PyVista mesh containing velocity field
        vel_field: Name of velocity field (default: 'U')

    Returns:
        Velocity magnitude array

    Raises:
        KeyError: If velocity field not found in mesh
    """
    if vel_field not in mesh.array_names:
        raise KeyError(
            f"Velocity field '{vel_field}' not found. Available: {mesh.array_names}")

    U = mesh[vel_field]
    return np.linalg.norm(U, axis=1) if U.ndim > 1 else U


def get_mesh_bounds(mesh: pv.UnstructuredGrid) -> dict:
    """
    Get mesh bounding box information.

    Args:
        mesh: PyVista mesh

    Returns:
        Dictionary with x, y, z min/max values
    """
    bounds = mesh.bounds
    return {
        'x_min': bounds[0],
        'x_max': bounds[1],
        'y_min': bounds[2],
        'y_max': bounds[3],
        'z_min': bounds[4],
        'z_max': bounds[5]
    }


def extract_surface_data(
    mesh: pv.UnstructuredGrid,
    field_name: str
) -> np.ndarray:
    """
    Extract field data from mesh surface.

    Args:
        mesh: PyVista mesh
        field_name: Name of field to extract (e.g., 'p', 'U', 'Cp')

    Returns:
        Extracted field data

    Raises:
        KeyError: If field not found in mesh
    """
    if field_name not in mesh.array_names:
        raise KeyError(f"Field '{field_name}' not found. Available: {mesh.array_names}")

    surface = mesh.extract_surface()
    return surface[field_name]


def find_latest_force_coeffs_file(case_dir: Path) -> Path:
    """
    Find the most recent OpenFOAM forceCoeffs output file for a given case.

    Args:
        case_dir (Path): Path to the case directory.

    Returns:
        Path: Path to the selected coefficient.dat file.

    Raises:
        FileNotFoundError: If the forceCoeffs directory or file cannot be found.
    """
    force_coeffs_dir = case_dir / "postProcessing" / "forceCoeffs"
    candidates = list(force_coeffs_dir.rglob("coefficient.dat"))
    if not candidates:
        raise FileNotFoundError(f"No coefficient.dat found under: {force_coeffs_dir}")

    times = []
    for p in candidates:
        t = float(p.parent.name)
        times.append((t, p.stat().st_mtime, p))
    times.sort(reverse=True)
    return times[0][2]


def read_force_coeffs_dat(file_path: Path) -> dict[str, np.ndarray]:
    """
    Read an OpenFOAM coefficient.dat file into named columns.

    Args:
        file_path (Path): Path to coefficient.dat.

    Returns:
        dict[str, np.ndarray]: Mapping from column name to data array.

    Raises:
        ValueError: If the file cannot be parsed into numeric data.
    """
    header_cols: Optional[List[str]] = None
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s.startswith("#"):
                continue

            maybe = s.lstrip("#").strip()
            parts = maybe.split()
            if parts and parts[0] == "Time":
                header_cols = parts

    data = np.loadtxt(file_path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if header_cols is None:
        header_cols = ["Time", "Cd", "Cs", "Cl", "CmRoll", "CmPitch", "CmYaw"]

    ncols = min(len(header_cols), data.shape[1])
    return {header_cols[i]: data[:, i] for i in range(ncols)}


def combine_postprocessing_summaries(
        working_path: Path,
        output_filename: str = "combined_postprocessing_summary.csv"
) -> None:
    """
    Combine all postprocessing_summary.csv files from subdirectories of working_path
    into a single CSV, sorted by AngleOfAttack.

    Args:
        working_path (Path): Path containing case directories with
            postprocessing_summary.csv files.
        output_filename (str): Name of the combined CSV file to create.
    """
    working_path = Path(working_path)
    summary_files = list(working_path.rglob("postprocessing_summary.csv"))
    if not summary_files:
        logger.log(f"No postprocessing_summary.csv files found in {working_path}")
        return

    dfs = []
    for file in summary_files:
        try:
            df = pd.read_csv(file)
            df["CaseDir"] = str(file.parent)
            dfs.append(df)
        except Exception as e:
            logger.log(f"Failed to read {file}: {e}")

    if not dfs:
        logger.log("No valid summary files to combine.")
        return

    combined_df = pd.concat(dfs, ignore_index=True)
    # Sort by AngleOfAttack if present
    if "AngleOfAttack" in combined_df.columns:
        combined_df = combined_df.sort_values("AngleOfAttack")
    output_path = working_path / output_filename
    combined_df.to_csv(output_path, index=False)
    logger.log(f"Combined summary saved to {output_path}")


def export_postprocessing_summary(
    case_dir: Path,
    output_csv: Path = None,
    mean_n: int = 1
) -> None:
    """
    Export postprocessing summary CSV for an OpenFOAM airfoil case.

    Args:
        case_dir (Path): Path to the case directory containing the results.
        output_csv (Path): Output CSV file path (default:
                <case_dir>/postprocessing_summary.csv).
        mean_n (int): Number of last timesteps to average for coefficients (default: 1).
    """
    import json
    import numpy as np

    output_csv = output_csv or (case_dir / "postprocessing_summary.csv")
    json_file = case_dir / f"{case_dir.name}.json"
    bc_csv = case_dir / "boundary_conditions_summary.csv"
    coeff_file = find_latest_force_coeffs_file(case_dir)
    coeffs = read_force_coeffs_dat(coeff_file)

    with open(json_file, "r") as f:
        j = json.load(f)
    designation = j["Airfoil"].get("Designation", "")
    aoa = j["Airfoil"].get("AngleOfAttack", "")

    df_bc = pd.read_csv(bc_csv)
    bc = dict(zip(df_bc['Parameter'], df_bc['Value']))
    mach = bc.get("Mach Number", "")
    reynolds = bc.get("Reynolds Number", "")
    velocity = bc.get("Velocity", "")
    density = bc.get("Density", "")
    chord = bc.get("Chord Length", "")

    n = len(coeffs["Time"])
    idx_start = max(0, n - mean_n)
    cd = np.mean(coeffs["Cd"][idx_start:])
    cl = np.mean(coeffs["Cl"][idx_start:])
    cmpitch = np.mean(coeffs.get("CmPitch", coeffs.get("Cm", np.zeros(n)))[idx_start:])

    df = pd.DataFrame([{
        "Designation": designation,
        "AngleOfAttack": aoa,
        "Cd": cd,
        "Cl": cl,
        "CmPitch": cmpitch,
        "Mach": mach,
        "Reynolds": reynolds,
        "Velocity": velocity,
        "Density": density,
        "Chord": chord
    }])

    df.to_csv(output_csv, index=False)
    logger.log(f"Exported summary to {output_csv}")


def read_uiuc_reference_csv(csv_path):
    """
    Reads a UIUC/Xfoil-style reference CSV, skipping header lines.
    Returns a pandas DataFrame with columns: Alpha, Cl, Cd, etc.
    """
    with open(csv_path, 'r') as f:
        lines = f.readlines()
    # Find the line where the data table starts
    for idx, line in enumerate(lines):
        if line.strip().startswith("Alpha,"):
            start_idx = idx
            break
    else:
        raise ValueError(f"Could not find data table in {csv_path}")
    # Read the data table into a DataFrame
    data_str = "".join(lines[start_idx:])
    df = pd.read_csv(StringIO(data_str))
    return df


def interpolate_airfoil_coefficients(
        csv_path: Path,
        angle_of_attack: float
) -> dict[str, float]:
    """
    Interpolate aerodynamic coefficients (Cl, Cd, Cm) from a UIUC/Xfoil-style CSV
    for a given angle of attack in the range -4 to +17 degrees.

    Args:
        csv_path (Path): Path to the reference CSV file.
        angle_of_attack (float): Angle of attack in degrees.

    Returns:
        dict[str, float]: Dictionary containing interpolated 'Cl', 'Cd', and 'Cm' values.

    Raises:
        ValueError: If angle of attack is outside the valid range (-4 to +17 degrees).
        FileNotFoundError: If CSV file not found.
        KeyError: If required columns are missing in CSV.
    """
    if not isinstance(csv_path, Path):
        csv_path = Path(csv_path)

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    if angle_of_attack < -4 or angle_of_attack > 17:
        raise ValueError(
            f"Angle of attack {angle_of_attack}° is outside valid range [-4°, +17°]"
        )

    # Read the CSV file
    df = read_uiuc_reference_csv(csv_path)

    # Ensure required columns exist
    required_cols = ['Alpha', 'Cl', 'Cd', 'Cm']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing columns in CSV: {missing_cols}")

    # Extract alpha values and coefficients
    alpha = df['Alpha'].values
    cl = df['Cl'].values
    cd = df['Cd'].values
    cm = df['Cm'].values

    # Create interpolation functions
    interp_cl = interp1d(alpha, cl, kind='cubic', fill_value='extrapolate')
    interp_cd = interp1d(alpha, cd, kind='cubic', fill_value='extrapolate')
    interp_cm = interp1d(alpha, cm, kind='cubic', fill_value='extrapolate')

    # Interpolate at the requested angle of attack
    cl_interp = float(interp_cl(angle_of_attack))
    cd_interp = float(interp_cd(angle_of_attack))
    cm_interp = float(interp_cm(angle_of_attack))

    return {
        'Cl': cl_interp,
        'Cd': cd_interp,
        'Cm': cm_interp,
        'Alpha': angle_of_attack
    }
