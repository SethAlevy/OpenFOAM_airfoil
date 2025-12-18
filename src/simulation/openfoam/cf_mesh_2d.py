from pathlib import Path
import numpy as np
from airfoil.airfoil import Airfoil
from templates.initial_settings_template import Settings
from simulation.openfoam.block_mesh import get_bounding_box
from utils.utilities import export_domain_to_fms, create_stl_bounding_box


def local_refinement(refinement_dict: dict) -> str:
    """
    Generate local refinement section for meshDict which is dedicated for patches.

    Args:
        refinement_dict (dict): Dictionary with patch names as keys and refinement
            settings as values. Each value should be a dict with level and thickness.

    Returns:
        str: Containing the localRefinement section.
    """
    if not refinement_dict:
        return ""

    refinement_entries = []
    for patch_name, settings in refinement_dict.items():
        level = settings.get("Level", 0)
        thickness = settings.get("Thickness", 0.05)

        entry = f"""    {patch_name}
    {{
        additionalRefinementLevels {level};
        refinementThickness {thickness};
    }}"""
        refinement_entries.append(entry)

    refinement_content = "\n\n".join(refinement_entries)

    return f"""
localRefinement
{{
{refinement_content}
}}
"""


def boundary_layers(layers_dict: dict) -> str:
    """
    Generate boundary layers section for meshDict.

    Args:
        layers_dict (dict): Dictionary with patch names as keys and layer settings
            as values. Each value should be a dict with: NLayers, ThicknessRatio,
            MaxFirstLayerThickness, AllowedDiscontinuity.

    Returns:
        str: Containing the boundaryLayers section.
    """
    if not layers_dict:
        return ""

    layer_entries = []
    for patch_name, settings in layers_dict.items():
        n_layers = settings.get("NLayers", 3)
        thickness_ratio = settings.get("ThicknessRatio", 1.2)
        max_first_layer = settings.get("MaxFirstLayerThickness", 0.001)
        discontinuity = settings.get("AllowedDiscontinuity", 1)

        entry = f"""    {patch_name}
    {{
        nLayers              {n_layers};
        thicknessRatio       {thickness_ratio};
        maxFirstLayerThickness {max_first_layer};
        allowedDiscontinuity {discontinuity};
    }}"""
        layer_entries.append(entry)

    layer_content = "\n\n".join(layer_entries)

    return f"""// Boundary layer settings
boundaryLayers
{{
    patchBoundaryLayers
    {{
    {layer_content}
    }}
}}
"""


def sphere_refinement(
    center: tuple[float, float, float],
    radius: float,
    cell_size: int,
    name: str
) -> str:
    """
    Generate sphere refinement section for meshDict.

    Args:
        center (tuple): Center of the sphere (x, y, z).
        radius (float): Radius of the sphere.
        cell_size (int): Desired cell size within the sphere.
        name (str): Name identifier for the refinement object.

    Returns:
        str: Containing the sphere refinement entry.
    """
    return f"""    {name}
    {{
        type            sphere;
        centre          ({center[0]} {center[1]} {center[2]});
        radius          {radius};
        cellSize        {cell_size};
    }}\n"""


def box_refinement(
    min_corner: tuple[float, float, float],
    max_corner: tuple[float, float, float],
    cell_size: int,
    name: str
) -> str:
    """
    Generate box refinement section for meshDict.

    Args:
        min_corner (tuple): Minimum corner of the box (x_min, y_min, z_min).
        max_corner (tuple): Maximum corner of the box (x_max, y_max, z_max).
        cell_size (int): Desired cell size within the box.
        name (str): Name identifier for the refinement object.

    Returns:
        str: Containing the box refinement entry.
    """
    box_center = (
        (min_corner[0] + max_corner[0]) / 2,
        (min_corner[1] + max_corner[1]) / 2,
        (min_corner[2] + max_corner[2]) / 2,
    )

    x_length = max_corner[0] - min_corner[0]
    y_length = max_corner[1] - min_corner[1]

    return f"""    {name}
    {{
        type            box;
        centre         ({box_center[0]} {box_center[1]} {box_center[2]});
        lengthX        {x_length};
        lengthY        {y_length};
        lengthZ        1;
        cellSize       {cell_size};
    }}\n"""


def object_refinements(airfoil: Airfoil, refinement_dict: dict) -> str:
    """
    Generate object refinements section for meshDict.

    Args:
        airfoil (Airfoil): Airfoil object containing airfoil properties.
        refinement_dict (dict): Dictionary with object names as keys and refinement
            settings as values. Each value should contain Type, and type-specific
            parameters like Center, Radius, CellSize for spheres or MinCorner,
            MaxCorner, CellSize for boxes.

    Returns:
        str: Containing the objectRefinements section.
    """
    if not refinement_dict:
        return ""
    
    objects_dicts = ""
    for key in refinement_dict.keys():
        if refinement_dict[key]["Type"] == "sphere":
            if "leadingEdge" in key.lower():
                centre = (0.0, 0.0, 0.0)
            elif "tip" in key.lower():
                centre = (
                    airfoil.upper_surface[0][-1],
                    airfoil.upper_surface[1][-1],
                    0.0
                )
            else:
                centre = refinement_dict[key].get("Center", [0.0, 0.0, 0.0])
            objects_dicts += sphere_refinement(
                center=centre,
                radius=refinement_dict[key].get("Radius", 0.2),
                cell_size=refinement_dict[key].get("CellSize", 0.01),
                name=key
            )
        elif refinement_dict[key]["Type"] == "box":
            x_min = refinement_dict[key].get("XMin", -0.5)
            x_max = refinement_dict[key].get("XMax", airfoil.chord * 2)
            y_min = refinement_dict[key].get("YMin", -0.5)
            y_max = refinement_dict[key].get("YMax", 0.5)

            objects_dicts += box_refinement(
                min_corner=refinement_dict[key].get(
                    "MinCorner",
                    [x_min, y_min, 0.0]
                ),
                max_corner=refinement_dict[key].get(
                    "MaxCorner",
                    [x_max, y_max, 0.0]
                ),
                cell_size=refinement_dict[key].get("CellSize", 0.01),
                name=key
            )

    return f"""
objectRefinements
{{
    // Object refinement settings can be added here
{objects_dicts}
}}
"""


def generate_cf_meshDict(
    airfoil: Airfoil,
    location_in_mesh: tuple[float, float],
    max_cell_size: float,
    local_refinement_dict: dict = None,
    boundary_layers_dict: dict = None,
    object_refinement_dict: dict = None
) -> str:
    """
    Generate meshDict for cartesian2DMesh with optional refinements and boundary layers.

    Args:
        airfoil (Airfoil): Airfoil object containing airfoil properties.
        location_in_mesh (tuple): Point inside the fluid domain (x, y).
        max_cell_size (float): Maximum cell size in the domain.
        local_refinement_dict (dict): Dictionary of patch names to refinement levels.
        boundary_layers_dict (dict): Dictionary of patch names to boundary layer
            settings.
        object_refinement_dict (dict): Dictionary of object names to refinement
            settings.

    Returns:
        str: Containing the complete meshDict content.
    """
    # Generate sections
    refinement_section = (
        local_refinement(local_refinement_dict)
        if local_refinement_dict else ""
    )
    layers_section = (
        boundary_layers(boundary_layers_dict)
        if boundary_layers_dict else ""
    )
    object_section = (
        object_refinements(airfoil, object_refinement_dict)
        if object_refinement_dict else ""
    )

    return f"""FoamFile
{{
    version   2.0;
    format    ascii;
    class     dictionary;
    location  "system";
    object    meshDict;
}}

surfaceFile         "constant/triSurface/domain.fms";
maxCellSize         {max_cell_size * 0.99999};
locationInMesh      ({location_in_mesh[0]} {location_in_mesh[1]});

{refinement_section}
{layers_section}
{object_section}
"""


def create_stl_domain(
    tri_surface_dir: Path,
    edge_patch_names: dict[str, str]
) -> None:
    """
    Create combined domain STL file from individual patch STL files.

    Args:
        tri_surface_dir (Path): Directory containing STL files.
        edge_patch_names (dict): Dictionary mapping edge types to patch names.

    Raises:
        FileNotFoundError: If airfoil.stl is not found in the directory.
    """
    airfoil_stl_path = tri_surface_dir / "airfoil.stl"
    if not airfoil_stl_path.exists():
        raise FileNotFoundError(f"Airfoil STL not found at {airfoil_stl_path}")

    stl_files_to_combine = [airfoil_stl_path] + [
        tri_surface_dir / f"{name}.stl" 
        for name in edge_patch_names.values()
    ]

    with open(tri_surface_dir / "domain.stl", "wb") as combined_file:
        for stl_file in stl_files_to_combine:
            if stl_file.exists():
                with open(stl_file, "rb") as individual_file:
                    combined_file.write(individual_file.read())


def cf_mesh_dict(
    airfoil: Airfoil,
    setup: Settings,
    system_dir: Path,
    stl_dir: Path = None
) -> None:
    """
    Generate .fms file and meshDict for the cartesian2DMesh workflow.

    Creates individual STL files for boundary patches and a combined domain.stl
    file for visualization. Also generates the .fms file describing the 2D domain
    geometry and the meshDict configuration file.

    Args:
        airfoil (Airfoil): Airfoil object containing geometry data.
        setup (Settings): Settings object containing mesh parameters.
        system_dir (Path): Path to the OpenFOAM case's system directory.
        stl_dir (Path): Path where the airfoil STL file is located and where to
            save generated STLs. If not given, defaults to case_dir/constant/triSurface.
    """
    if stl_dir is None:
        case_dir = system_dir.parent
        stl_dir = case_dir / "constant" / "triSurface"

    cf_mesh_settings = setup.mesh_settings.get("CfMesh", {})

    x_min, x_max, y_min, y_max, z_min, z_max = get_bounding_box(
        airfoil, setup.mesh_settings.get("BoundingBox", {})
    )

    boundary_conditions = setup.simulation_settings.get("BoundaryConditions", {})

    edge_patch_names = {
        "inlet": boundary_conditions.get("InletPatchName", "inlet"),
        "outlet": boundary_conditions.get("OutletPatchName", "outlet"),
        "lowerWall": boundary_conditions.get("LowerWallPatchName", "lowerWall"),
        "upperWall": boundary_conditions.get("UpperWallPatchName", "upperWall"),
    }

    create_stl_bounding_box(
        x_min, x_max, y_min, y_max, z_min, z_max,
        edge_patch_names, stl_dir
    )

    create_stl_domain(stl_dir, edge_patch_names)

    x_contour = np.concatenate(
        (airfoil.upper_surface[0], airfoil.lower_surface[0][::-1])
    )
    y_contour = np.concatenate(
        (airfoil.upper_surface[1], airfoil.lower_surface[1][::-1])
    )

    export_domain_to_fms(
        x_contour, y_contour,
        x_min, x_max, y_min, y_max, z_min, z_max,
        patch_names=edge_patch_names,
        output_path=stl_dir / "domain.fms"
    )

    location_in_mesh = (x_min + 0.1 * airfoil.chord, 0.0)

    mesh_dict_content = generate_cf_meshDict(
        airfoil=airfoil,
        location_in_mesh=location_in_mesh,
        max_cell_size=cf_mesh_settings.get("MaxCellSize", 0.02),
        local_refinement_dict=cf_mesh_settings.get("LocalRefinement", {}),
        boundary_layers_dict=cf_mesh_settings.get("BoundaryLayers", {}),
        object_refinement_dict=cf_mesh_settings.get("ObjectRefinements", {})
    )
    with open(system_dir / "meshDict", "w") as f:
        f.write(mesh_dict_content)
