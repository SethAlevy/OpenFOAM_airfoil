"""
Utilities for geometric manipulations and mesh exports.
"""

import numpy as np
import pyvista as pv
from pathlib import Path
from stl import mesh as np_mesh, Mode
from scipy.interpolate import interp1d
from utils.logger import SimpleLogger as logger


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


def rotate_by_alpha(alpha: float,
                    x: np.ndarray,
                    y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Rotate coordinates by angle of attack alpha.

    Args:
        alpha (float): angle of attack in degrees.
        x (np.ndarray): x-coordinates.
        y (np.ndarray): y-coordinates.

    Returns:
        tuple[np.ndarray, np.ndarray]: rotated (x_rotated, y_rotated) coordinates.
    """
    alpha_rad = np.radians(alpha)
    x_rot = x * np.cos(alpha_rad) - y * np.sin(alpha_rad)
    y_rot = x * np.sin(alpha_rad) + y * np.cos(alpha_rad)
    return x_rot, y_rot


def resample_line(x: np.ndarray,
                  y: np.ndarray,
                  n_points: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Resample a line defined by x and y coordinates to have n_points.

    Args:
        x (np.ndarray): Original x-coordinates.
        y (np.ndarray): Original y-coordinates.
        n_points (int): Number of points to resample to.

    Returns:
        tuple[np.ndarray, np.ndarray]: Resampled x and y coordinates.
    """
    interp_func = interp1d(x, y, kind="linear")
    x_new = np.linspace(np.min(x), np.max(x), n_points)
    y_new = interp_func(x_new)
    return x_new, y_new


def resample_by_arc_length(x: np.ndarray,
                           y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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
    dx = np.diff(x)
    dy = np.diff(y)
    segment_lengths = np.sqrt(dx**2 + dy**2)
    cumulative_length = np.concatenate(([0], np.cumsum(segment_lengths)))

    interp_x = interp1d(cumulative_length, x, kind="linear")
    interp_y = interp1d(cumulative_length, y, kind="linear")

    total_length = cumulative_length[-1]
    uniform_arc_lengths = np.linspace(0, total_length, n_points)

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
    resample: bool = True,
) -> None:
    """
    Creates and saves a 3D extruded closed airfoil as a numpy-stl ASCII file.
    Optionally resamples coordinates to have constant arc length spacing for
    better mesh quality.

    Args:
        x_upper (np.ndarray): x-coordinates of the upper surface (LE to TE).
        y_upper (np.ndarray): y-coordinates of the upper surface.
        x_lower (np.ndarray): x-coordinates of the lower surface (TE to LE).
        y_lower (np.ndarray): y-coordinates of the lower surface.
        output_path (Path): STL file path.
        thickness (float): Extrusion thickness in meters (default: 1e-3).
        resample (bool): If True, resample with uniform arc length spacing.
    """
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
    faces.extend([front_pts[0], front_pts[i + 1], front_pts[i + 2]]
                 for i in range(n - 2))
    faces.extend([back_pts[0], back_pts[i + 2], back_pts[i + 1]]
                 for i in range(n - 2))
    for i in range(n - 1):
        faces.extend(
            (
                [front_pts[i], front_pts[i + 1], back_pts[i + 1]],
                [front_pts[i], back_pts[i + 1], back_pts[i]],
            )
        )

    faces = np.array(faces)
    airfoil_mesh = np_mesh.Mesh(
        np.zeros(faces.shape[0], dtype=np_mesh.Mesh.dtype), name=b"airfoil"
    )
    airfoil_mesh.vectors = faces
    airfoil_mesh.save(str(output_path), mode=Mode.ASCII)


def create_stl_bounding_box(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    z_min: float,
    z_max: float,
    patch_names: dict,
    output_dir: Path,
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
        patch_names (dict): Mapping patch roles to names.
        output_dir (Path): Directory to save the STL files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    def make_face(v1, v2, v3, v4, name):
        face_mesh = np_mesh.Mesh(
            np.zeros(2, dtype=np_mesh.Mesh.dtype), name=name.encode()
        )
        face_mesh.vectors[0] = [v1, v2, v3]
        face_mesh.vectors[1] = [v1, v3, v4]
        face_mesh.save(str(output_dir / f"{name}.stl"), mode=Mode.ASCII)

    p0 = (x_min, y_min, z_min)
    p1 = (x_max, y_min, z_min)
    p2 = (x_max, y_max, z_min)
    p3 = (x_min, y_max, z_min)
    p4 = (x_min, y_min, z_max)
    p5 = (x_max, y_min, z_max)
    p6 = (x_max, y_max, z_max)
    p7 = (x_min, y_max, z_max)

    make_face(p3, p7, p4, p0, patch_names["inlet"])
    make_face(p1, p5, p6, p2, patch_names["outlet"])
    make_face(p0, p4, p5, p1, patch_names["lowerWall"])
    make_face(p2, p6, p7, p3, patch_names["upperWall"])


def export_airfoil_to_stl_2d(
    x_upper: np.ndarray,
    y_upper: np.ndarray,
    x_lower: np.ndarray,
    y_lower: np.ndarray,
    output_path: Path,
    z_value: float = 0.0,
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

    if not (
        np.isclose(x_contour[0], x_contour[-1])
        and np.isclose(y_contour[0], y_contour[-1])
    ):
        x_contour = np.append(x_contour, x_contour[0])
        y_contour = np.append(y_contour, y_contour[0])
        n += 1

    with open(output_path, "w") as stl:
        stl.write("solid airfoil2d\n")
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
    output_path: Path,
) -> None:
    """
    Export the domain around the airfoil as an ASCII FMS file. This allows cfMesh
    to read the patch names for defined faces.

    Args:
        airfoil_x (np.ndarray): x-coordinates of the closed airfoil contour.
        airfoil_y (np.ndarray): y-coordinates of the closed airfoil contour.
        x_min (float): Minimum x-coordinate of the bounding box.
        x_max (float): Maximum x-coordinate of the bounding box.
        y_min (float): Minimum y-coordinate of the bounding box.
        y_max (float): Maximum y-coordinate of the bounding box.
        z_front (float): z-coordinate of the front face.
        z_back (float): z-coordinate of the back face.
        patch_names (dict): Mapping patch roles to their names.
        output_path (Path): Path to save the FMS file.
    """
    all_patches_with_types = {"airfoil": "wall"}
    all_patches_with_types.update({name: "patch" for name in patch_names.values()})

    num_patches = len(all_patches_with_types)
    patch_to_index = {name: idx for idx, name in enumerate(all_patches_with_types)}

    vertices = []
    num_airfoil_verts = len(airfoil_x)
    vertices.extend((airfoil_x[i], airfoil_y[i], z_front)
                    for i in range(num_airfoil_verts))
    vertices.extend((airfoil_x[i], airfoil_y[i], z_back)
                    for i in range(num_airfoil_verts))

    bbox_front_start = 2 * num_airfoil_verts
    vertices.extend(
        [
            (x_min, y_max, z_front),
            (x_max, y_max, z_front),
            (x_max, y_min, z_front),
            (x_min, y_min, z_front),
        ]
    )

    bbox_back_start = bbox_front_start + 4
    vertices.extend(
        [
            (x_min, y_max, z_back),
            (x_max, y_max, z_back),
            (x_max, y_min, z_back),
            (x_min, y_min, z_back),
        ]
    )

    faces = []
    faces.extend(
        ((0, i, i + 1), patch_to_index["airfoil"])
        for i in range(1, num_airfoil_verts - 1)
    )
    for i in range(1, num_airfoil_verts - 1):
        v0 = num_airfoil_verts
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

    v0, v1, v2, v3 = (
        bbox_front_start + 0,
        bbox_front_start + 1,
        bbox_back_start + 1,
        bbox_back_start + 0,
    )
    faces.append(((v0, v1, v2), patch_to_index[patch_names["upperWall"]]))
    faces.append(((v0, v2, v3), patch_to_index[patch_names["upperWall"]]))

    v0, v1, v2, v3 = (
        bbox_front_start + 1,
        bbox_front_start + 2,
        bbox_back_start + 2,
        bbox_back_start + 1,
    )
    faces.append(((v0, v1, v2), patch_to_index[patch_names["outlet"]]))
    faces.append(((v0, v2, v3), patch_to_index[patch_names["outlet"]]))

    v0, v1, v2, v3 = (
        bbox_front_start + 2,
        bbox_front_start + 3,
        bbox_back_start + 3,
        bbox_back_start + 2,
    )
    faces.append(((v0, v1, v2), patch_to_index[patch_names["lowerWall"]]))
    faces.append(((v0, v2, v3), patch_to_index[patch_names["lowerWall"]]))

    v0, v1, v2, v3 = (
        bbox_front_start + 3,
        bbox_front_start + 0,
        bbox_back_start + 0,
        bbox_back_start + 3,
    )
    faces.append(((v0, v1, v2), patch_to_index[patch_names["inlet"]]))
    faces.append(((v0, v2, v3), patch_to_index[patch_names["inlet"]]))

    with open(output_path, "w") as f:
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
            f.write(f"(({face_verts[0]} {face_verts[1]} {face_verts[2]}) {patch_idx})\n")
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
