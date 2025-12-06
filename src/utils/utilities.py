import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
from stl import mesh as np_mesh, Mode  # Import Mode from the top-level stl module


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


def export_airfoil_to_stl_3d(
        x: np.ndarray,
        y: np.ndarray,
        output_path: Path,
        thickness: float = 1e-3
) -> None:
    """
    (Restored)
    Creates and saves a 3D extruded airfoil as a numpy-stl ASCII file.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    assert len(x) == len(y) and len(x) >= 3

    if not (np.isclose(x[0], x[-1]) and np.isclose(y[0], y[-1])):
        x = np.r_[x, x[0]]
        y = np.r_[y, y[0]]

    z_top = +0.5 * thickness
    z_bot = -0.5 * thickness
    n_pts = len(x)

    top_pts = np.c_[x, y, np.full(n_pts, z_top)]
    bot_pts = np.c_[x, y, np.full(n_pts, z_bot)]

    faces = []
    # Top cap
    for i in range(1, n_pts - 2):
        faces.append([top_pts[0], top_pts[i], top_pts[i + 1]])
    # Bottom cap
    for i in range(1, n_pts - 2):
        faces.append([bot_pts[0], bot_pts[i + 1], bot_pts[i]])
    # Side walls
    for i in range(n_pts - 1):
        p0_top, p1_top = top_pts[i], top_pts[i + 1]
        p0_bot, p1_bot = bot_pts[i], bot_pts[i + 1]
        faces.append([p0_top, p1_top, p1_bot])
        faces.append([p0_top, p1_bot, p0_bot])

    faces = np.array(faces)
    airfoil_mesh = np_mesh.Mesh(
        np.zeros(faces.shape[0], dtype=np_mesh.Mesh.dtype), name=b'airfoil')
    airfoil_mesh.vectors = faces

    # Save the mesh to the specified file path using the correct Mode enum
    airfoil_mesh.save(str(output_path), mode=Mode.ASCII)


def create_stl_bounding_box(
    x_min: float, x_max: float, y_min: float, y_max: float, z_min: float, z_max: float,
    patch_names: dict, output_dir: Path
) -> None:
    """
    (Modified)
    Creates and saves 6 STL files representing the faces of a bounding box.
    Each STL is named according to the patch name.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    def make_face(v1, v2, v3, v4, name):
        # The name of the solid inside the STL must match the patch name
        face_mesh = np_mesh.Mesh(
            np.zeros(2, dtype=np_mesh.Mesh.dtype), name=name.encode())
        face_mesh.vectors[0] = [v1, v2, v3]
        face_mesh.vectors[1] = [v1, v3, v4]
        # Save to a file with the same name
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
        x: np.ndarray,
        y: np.ndarray,
        output_path: Path,
        z_value: float = 0.0
) -> None:
    """
    Export a 2D airfoil contour as a flat ASCII STL file at a single z-value.

    Args:
        x (np.ndarray): x-coordinates of the closed airfoil contour.
        y (np.ndarray): y-coordinates of the closed airfoil contour.
        output_path (Path): STL file path.
        z_value (float): z-coordinate for all vertices (default: 0.0).
    """
    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)

    with open(output_path, 'w') as stl:
        stl.write("solid airfoil2d\n")
        for i in range(n - 1):
            v0 = (x[i], y[i], z_value)
            v1 = (x[i + 1], y[i + 1], z_value)
            v2 = (x[0], y[0], z_value)  # Reference vertex for triangle fan
            stl.write("  facet normal 0 0 1\n    outer loop\n")
            stl.write(f"      vertex {v0[0]} {v0[1]} {v0[2]}\n")
            stl.write(f"      vertex {v1[0]} {v1[1]} {v1[2]}\n")
            stl.write(f"      vertex {v2[0]} {v2[1]} {v2[2]}\n")
            stl.write("    endloop\n  endfacet\n")
        stl.write("endsolid airfoil2d\n")


def export_domain_to_fms(
    airfoil_x: np.ndarray, airfoil_y: np.ndarray,
    x_min: float, x_max: float, y_min: float, y_max: float,
    patch_names: dict, output_path: Path
) -> None:
    """
    (Corrected)
    Exports a 3D domain to .fms file format with properly separated
    airfoil and bounding box surfaces.
    """
    # Create the header section with correct types
    all_patches_with_types = {"airfoil": "wall"}
    all_patches_with_types.update({name: "patch" for name in patch_names.values()})

    num_patches = len(all_patches_with_types)
    patch_to_index = {name: idx for idx,
                      name in enumerate(all_patches_with_types.keys())}

    # Z coordinates for 3D extrusion
    z_front = 0.01
    z_back = -0.01

    # Build vertex list
    vertices = []

    # Airfoil vertices - front surface (z = z_front)
    num_airfoil_verts = len(airfoil_x)
    for i in range(num_airfoil_verts):
        vertices.append((airfoil_x[i], airfoil_y[i], z_front))

    # Airfoil vertices - back surface (z = z_back)
    for i in range(num_airfoil_verts):
        vertices.append((airfoil_x[i], airfoil_y[i], z_back))

    # Bounding box vertices - front
    bbox_front_start = 2 * num_airfoil_verts
    vertices.extend([
        (x_min, y_max, z_front),  # 0: Top-left front
        (x_max, y_max, z_front),  # 1: Top-right front
        (x_max, y_min, z_front),  # 2: Bottom-right front
        (x_min, y_min, z_front),  # 3: Bottom-left front
    ])

    # Bounding box vertices - back
    bbox_back_start = bbox_front_start + 4
    vertices.extend([
        (x_min, y_max, z_back),   # 0: Top-left back
        (x_max, y_max, z_back),   # 1: Top-right back
        (x_max, y_min, z_back),   # 2: Bottom-right back
        (x_min, y_min, z_back),   # 3: Bottom-left back
    ])

    # Build face list
    faces = []

    # 1. Airfoil front cap (fan triangulation from first vertex as centroid)
    for i in range(1, num_airfoil_verts - 1):
        faces.append(((0, i, i + 1), patch_to_index["airfoil"]))

    # 2. Airfoil back cap (reverse winding for correct normal)
    for i in range(1, num_airfoil_verts - 1):
        v0 = num_airfoil_verts  # First vertex on back
        v1 = num_airfoil_verts + i + 1
        v2 = num_airfoil_verts + i
        faces.append(((v0, v1, v2), patch_to_index["airfoil"]))

    # 3. Airfoil side walls (connecting front and back, forming the extruded thickness)
    for i in range(num_airfoil_verts):
        v0_front = i
        v1_front = (i + 1) % num_airfoil_verts
        v0_back = num_airfoil_verts + i
        v1_back = num_airfoil_verts + (i + 1) % num_airfoil_verts

        # Two triangles for each quad section
        faces.append(((v0_front, v1_front, v1_back), patch_to_index["airfoil"]))
        faces.append(((v0_front, v1_back, v0_back), patch_to_index["airfoil"]))

    # 4. Bounding box faces (each face is 2 triangles)
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

    # Write to file
    with open(output_path, 'w') as f:
        # --- 1. Write Boundary Patch Definitions ---
        f.write(f"{num_patches}\n(\n\n")
        for name, patch_type in all_patches_with_types.items():
            f.write(f"{name} {patch_type}\n\n")
        f.write(")\n\n\n")

        # --- 2. Write Vertices ---
        f.write(f"{len(vertices)}\n(\n")
        for v in vertices:
            f.write(f"({v[0]} {v[1]} {v[2]})\n")
        f.write(")\n\n\n")

        # --- 3. Write Faces ---
        f.write(f"{len(faces)}\n(\n")
        for face_verts, patch_idx in faces:
            f.write(
                f"(({face_verts[0]} {face_verts[1]} {face_verts[2]}) {patch_idx})\n")
        f.write(")\n\n\n")

        # --- 4. Write Edges (derived from faces) ---
        edges = set()
        for face_verts, _ in faces:
            for i in range(3):
                v1, v2 = face_verts[i], face_verts[(i + 1) % 3]
                edges.add((min(v1, v2), max(v1, v2)))

        f.write(f"{len(edges)}\n(\n")
        for e in sorted(edges):
            f.write(f"({e[0]} {e[1]})\n")
        f.write(")\n\n")

        # --- 5. Three empty sections ---
        f.write("0()\n")
        f.write("0()\n")
        f.write("0()\n")
