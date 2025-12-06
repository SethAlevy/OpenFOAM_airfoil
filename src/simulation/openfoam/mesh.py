from templates.airfoil_template import Airfoil
from templates.initial_settings_template import Settings
from pathlib import Path
import numpy as np
import subprocess
from stl import mesh as np_mesh, Mode  # Import Mode from the top-level stl module
from utils.utilities import create_stl_bounding_box
from utils.utilities import export_airfoil_to_stl_3d
from utils.utilities import export_domain_to_fms


def block_mesh_dict(airfoil: Airfoil, setup: Settings, output_path: Path) -> None:
    """
    Fill the blockMeshDict templates according to the selected simulation setup and
    save the file in the given path.

    Args:
        airfoil (Airfoil): Object with the generated airfoil geometry.
        setup (Settings): Object with the simulation settings.
        output_path (Path): Path where the blockMeshDict file will be saved.
    """
    bounding_box = setup.mesh_settings.get("BoundingBox", {})
    x_min, x_max, y_min, y_max = get_bounding_box(airfoil, bounding_box)
    cell_size = bounding_box["BaseCellSize"]

    boundary_conditions = setup.simulation_settings.get("BoundaryConditions", {})
    inlet_patch = boundary_conditions.get("InletPatchName", "inlet")
    outlet_patch = boundary_conditions.get("OutletPatchName", "outlet")
    lower_wall_patch = boundary_conditions.get("LowerWallPatchName", "lowerWall")
    upper_wall_patch = boundary_conditions.get("UpperWallPatchName", "upperWall")
    front_patch = boundary_conditions.get("FrontPatchName", "front")
    back_patch = boundary_conditions.get("BackPatchName", "back")

    content = generate_block_mesh_dict(
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        nx=int((x_max - x_min) / cell_size),
        ny=int((y_max - y_min) / cell_size),
        inlet_patch=inlet_patch,
        outlet_patch=outlet_patch,
        lower_wall_patch=lower_wall_patch,
        upper_wall_patch=upper_wall_patch,
        front_patch=front_patch,
        back_patch=back_patch
    )

    output_path.mkdir(parents=True, exist_ok=True)
    with open(output_path / "blockMeshDict", 'w') as file:
        file.write(content)


def get_bounding_box(
        airfoil: Airfoil,
        bounding_box: dict
) -> tuple[float, float, float, float]:
    """
    Calculate the bounding box of the airfoil geometry assuming some distances
    from the airfoil to the domain boundaries.

    Args:
        airfoil (Airfoil): Object with the generated airfoil geometry.
        bounding_box (dict): Dictionary with the bounding box distances.

    Returns:
        tuple[float, float, float, float]: x_min, x_max, y_min, y_max of the
            bounding box.
    """
    chord = airfoil.chord

    # surface coords are given as 2D arrays (x_coords, y_coords)
    x_min = min(airfoil.upper_surface[0]) - bounding_box["InletDistance"] * chord
    x_max = max(airfoil.upper_surface[0]) + bounding_box["OutletDistance"] * chord
    y_min = min(airfoil.lower_surface[1]) - bounding_box["BottomDistance"] * chord
    y_max = max(airfoil.upper_surface[1]) + bounding_box["TopDistance"] * chord
    return x_min, x_max, y_min, y_max


def generate_block_mesh_dict(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    nx: int,
    ny: int,
    z_min: float = -0.0015,
    z_max: float = 0.0015,
    nz: int = 1,
    inlet_patch: str = "inlet",
    outlet_patch: str = "outlet",
    lower_wall_patch: str = "lowerWall",
    upper_wall_patch: str = "upperWall",
    front_patch: str = "front",
    back_patch: str = "back"
) -> str:
    """
    Generate a simple BlockMesh dict for a 2D airfoil case with customizable patch
    names.

    Args:
        x_min (float): Minimum x coordinate of the bounding box.
        x_max (float): Maximum x coordinate of the bounding box.
        y_min (float): Minimum y coordinate of the bounding box.
        y_max (float): Maximum y coordinate of the bounding box.
        nx (int): Number of cells in the x direction.
        ny (int): Number of cells in the y direction.
        z_min (float): Minimum z coordinate of the bounding box.
        z_max (float): Maximum z coordinate of the bounding box.
        nz (int): Number of cells in the z direction.
        inlet_patch (str): Name of the inlet patch.
        outlet_patch (str): Name of the outlet patch.
        lower_wall_patch (str): Name of the lower wall patch.
        upper_wall_patch (str): Name of the upper wall patch.
        front_patch (str): Name of the front patch.
        back_patch (str): Name of the back patch.

    Returns:
        str: The content of the blockMeshDict file.
    """
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

convertToMeters 1.0;

vertices
(
    ({x_min} {y_min} {z_min})
    ({x_max} {y_min} {z_min})
    ({x_max} {y_max} {z_min})
    ({x_min} {y_max} {z_min})
    ({x_min} {y_min} {z_max})
    ({x_max} {y_min} {z_max})
    ({x_max} {y_max} {z_max})
    ({x_min} {y_max} {z_max})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    {inlet_patch}
    {{
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }}
    {outlet_patch}
    {{
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }}
    {lower_wall_patch}
    {{
        type patch;
        faces
        (
            (0 1 5 4)
        );
    }}
    {upper_wall_patch}
    {{
        type patch;
        faces
        (
            (3 7 6 2)
        );
    }}
    {front_patch}
    {{
        type symmetryPlane;
        faces
        (
            (0 3 2 1)
        );
    }}
    {back_patch}
    {{
        type symmetryPlane;
        faces
        (
            (4 5 6 7)
        );
    }}
);

mergePatchPairs
(
);
"""
    return content


def snappy_hex_mesh_dict(airfoil: Airfoil, setup: Settings, output_path: Path) -> None:
    """
    Fill the snappyHexMeshDict templates according to the selected simulation setup and
    save the file in the given path.

    Args:
        airfoil (Airfoil): Object with the generated airfoil geometry.
        setup (Settings): Object with the simulation settings.
        output_path (Path): Path where the snappyHexMeshDict file will be saved.
    """
    snappy_hex_mesh_settings = setup.mesh_settings.get("SnappyHexMesh", {})

    # Geometry refinement regions
    refinement_box_level = snappy_hex_mesh_settings.get(
        "Geometry", {}).get("RefinementBoxLevel", 0)
    sphere_tip_level = snappy_hex_mesh_settings.get(
        "Geometry", {}).get("SphereTipLevel", 0)
    sphere_tip_radius_scale = snappy_hex_mesh_settings.get(
        "Geometry", {}).get("SphereTipRadius", 0)
    sphere_leading_level = snappy_hex_mesh_settings.get(
        "Geometry", {}).get("SphereLeadingEdgeLevel", 0)
    sphere_leading_radius_scale = snappy_hex_mesh_settings.get(
        "Geometry", {}).get("SphereLeadingEdgeRadius", 0)
    if refinement_box_level > 0:
        refinement_box = refinement_box_dict(airfoil)
    else:
        refinement_box = ''
    if sphere_tip_level > 0:
        sphere_tip = refine_sphere_dict(airfoil, "tip", sphere_tip_radius_scale)
    else:
        sphere_tip = ''
    if sphere_leading_level > 0:
        sphere_leading = refine_sphere_dict(
            airfoil, "leading", sphere_leading_radius_scale)
    else:
        sphere_leading = ''

    # Section dictionaries from JSON
    add_layers_controls = snappy_hex_mesh_settings.get("AddLayersControls", {})
    snap_controls = snappy_hex_mesh_settings.get("SnapControls", {})
    castellated_controls = snappy_hex_mesh_settings.get("CastellatedMeshControls", {})

    # Feature refinement level
    feature_refinement_level = castellated_controls.get("FeatureRefinementLevel", 0)
    content = generate_snappy_hex_mesh_dict(
        airfoil=airfoil,
        refinement_box=refinement_box,
        refinement_box_level=refinement_box_level,
        sphere_tip=sphere_tip,
        sphere_tip_level=sphere_tip_level,
        sphere_leading=sphere_leading,
        sphere_leading_level=sphere_leading_level,
        feature_refinement_level=feature_refinement_level,
        layer_controls=add_layers_controls,
        snap_controls=snap_controls,
        castellated_controls=castellated_controls,
    )

    output_path.mkdir(parents=True, exist_ok=True)
    with open(output_path / "snappyHexMeshDict", 'w') as file:
        file.write(content)


def refinement_box_dict(Airfoil: Airfoil) -> dict:
    """
    Define the refinement box around the airfoil. Size is assumed by taking the
    airfoils extreme coordinates and adding one chord length distance in each direction.

    Args:
        Airfoil (Airfoil): Object with the generated airfoil geometry.

    Returns:
        dict: Dictionary with the minimum and maximum coordinates of the refinement box.
    """
    chord = Airfoil.chord
    upper_surface = Airfoil.upper_surface
    lower_surface = Airfoil.lower_surface

    return {
        "min": [min(lower_surface[0]) - chord, min(lower_surface[1]) - chord, -0.01],
        "max": [max(upper_surface[0]) + chord, max(upper_surface[1]) + chord, 0.01],
    }


def refine_sphere_dict(Airfoil: Airfoil, position: str, radius_scale: float) -> dict:
    """
    Define a refinement sphere around a specific position of the airfoil.

    Args:
        Airfoil (Airfoil): Object with the generated airfoil geometry.
        position (str): Position on the airfoil to center the sphere (available
            "tip" or "leading").
        radius_scale (float): Scale factor for the sphere radius based on chord length.
    Returns:
        dict: Dictionary with the center coordinates and radius of the refinement
            sphere.
    """
    upper_surface = Airfoil.upper_surface
    chord = Airfoil.chord

    if position == "tip":
        center = [upper_surface[0][-1], upper_surface[1][-1], 0.0]
    elif position == "leading":
        center = [upper_surface[0][0], upper_surface[1][0], 0.0]

    return {
        "center": center,
        "radius": radius_scale * chord
    }


def get_location_in_mesh(airfoil: Airfoil) -> np.ndarray:
    """
    Get the location of the airfoil in the mesh for snappyHexMesh refinement from
    mean camber line.

    Args:
        airfoil (Airfoil): Object with the generated airfoil geometry.

    Returns:
        np.ndarray: Array with the x and y coordinates of the airfoil location.
    """

    x_coords = airfoil.mean_camber_line[0]
    y_coords = airfoil.mean_camber_line[1]
    locations = np.column_stack((x_coords, y_coords))
    return locations[10]


def generate_snappy_hex_mesh_dict(
    airfoil: Airfoil,
    refinement_box: dict = '',
    refinement_box_level: int = 0,
    sphere_tip: dict = '',
    sphere_tip_level: int = 0,
    sphere_leading: dict = '',
    sphere_leading_level: int = 0,
    feature_refinement_level: int = 0,
    layer_controls: dict = {},
    snap_controls: dict = {},
    castellated_controls: dict = {},
    mesh_quality_controls: dict = {},
):
    """
    Generate a snappyHexMeshDict for a semi-2D airfoil case, using helper functions
    for sections. Allows full control over the meshing process through settings
    file.

    Args:
        airfoil (Airfoil): Object with the generated airfoil geometry.
        refinement_box (dict): Dictionary with the minimum and maximum coordinates of
            the airfoil refinement box.
        refinement_box_level (int): Refinement level for the airfoil refinement box.
        sphere_tip (dict): Dictionary with the center coordinates and radius of the
            tip refinement sphere.
        sphere_tip_level (int): Refinement level for the tip refinement sphere.
        sphere_leading (dict): Dictionary with the center coordinates and radius of the
            leading edge refinement sphere.
        sphere_leading_level (int): Refinement level for the leading edge refinement
            sphere.
        feature_refinement_level (int): Refinement level for the edge features.
        layer_controls (dict): Dictionary with the addLayersControls section settings.
        snap_controls (dict): Dictionary with the snapControls section settings.
        castellated_controls (dict): Dictionary with the castellatedMeshControls
            section settings.
        mesh_quality_controls (dict): Dictionary with the meshQualityControls section
            settings.

    Returns:
        str: The content of the snappyHexMeshDict file.
    """
    ref_box_geo, ref_box_reg = box_geometry_and_refinement_str(
        "refinementBox", refinement_box, refinement_box_level)
    tip_sphere_geo, tip_sphere_reg = sphere_geometry_and_refinement_str(
        "tipSphere", sphere_tip, sphere_tip_level)
    lead_sphere_geo, lead_sphere_reg = sphere_geometry_and_refinement_str(
        "leadingEdgeSphere", sphere_leading, sphere_leading_level)

    features_str = feature_str(feature_refinement_level)

    geometry_section = f"""
    airfoil.stl
    {{
        type triSurfaceMesh;
        name airfoil;
    }}
    {ref_box_geo}
    {tip_sphere_geo}
    {lead_sphere_geo}
"""
    chord = airfoil.chord

    refinement_regions_str = "    refinementRegions\n    {\n"
    refinement_regions_str += ref_box_reg
    refinement_regions_str += tip_sphere_reg
    refinement_regions_str += lead_sphere_reg
    refinement_regions_str += "    }\n"

    castellated_mesh_controls_section = castellated_mesh_controls_str(
        chord,
        castellated_controls,
        features_str,
        refinement_regions_str
    )

    snap_controls_section = snap_controls_str(snap_controls)

    layer_controls_section = add_layers_controls_str(
        layer_controls,
    )

    mesh_quality_controls_section = mesh_quality_controls_str(mesh_quality_controls)

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}

castellatedMesh true;
snap            true;
addLayers       {'true' if layer_controls and layer_controls.get('NSurfaceLayers', 0) > 0 else 'false'};

geometry
{{
{geometry_section}
}}

{castellated_mesh_controls_section}
{snap_controls_section}
{layer_controls_section}
{mesh_quality_controls_section}
debug 0;
mergeTolerance 1E-6;
"""
    return content


def box_geometry_and_refinement_str(
        name: str,
        box: dict,
        level: int
) -> tuple[str, str]:
    """
    Generate geometry and refinement region strings for a refinement box around the
    airfoil.

    Args:
        name (str): Name of the refinement box.
        box (dict): Dictionary with the minimum and maximum coordinates of the box.
        level (int): Refinement level for the box.

    Returns:
        tuple[str, str]: Geometry string and refinement region string for the box.
    """
    if not box:
        return "", ""
    geometry_str = f"""
    {name}
    {{
        type searchableBox;
        min ({box['min'][0]} {box['min'][1]} {box['min'][2]});
        max ({box['max'][0]} {box['max'][1]} {box['max'][2]});
    }}
    """
    refinement_str = f"""        {name}
        {{
            mode inside;
            levels ((1E15 {level}));
        }}
"""
    return geometry_str, refinement_str


def sphere_geometry_and_refinement_str(
        name: str,
        sphere: dict,
        level: int
) -> tuple[str, str]:
    """
    Generate geometry and refinement region strings for a sphere (around the tip or
    leading edge).

    Args:
        name (str): Name of the refinement sphere.
        sphere (dict): Dictionary with the center coordinates and radius of the sphere.
        level (int): Refinement level for the sphere.

    Returns:
        tuple[str, str]: Geometry string and refinement region string for the sphere.
    """
    if not sphere:
        return "", ""
    geometry_str = f"""
    {name}
    {{
        type searchableSphere;
        centre ({sphere['center'][0]} {sphere['center'][1]} {sphere['center'][2]});
        radius {sphere['radius']};
    }}
    """
    refinement_str = f"""        {name}
        {{
            mode inside;
            levels ((1E15 {level}));
        }}
"""
    return geometry_str, refinement_str


def feature_str(feature_level: int = 0, feature_file: str = "airfoil.eMesh"):
    """
    Generate the features section string for snappyHexMeshDict. If level is 0 skip
    the section.

    Args:
        feature_level (int): Refinement level for the edge features.
        feature_file (str): Filename of the edge feature file.
    """
    if feature_level != 0:
        return f"""    features
    (
        {{
            file "{feature_file}";
            level {feature_level};
        }}
    );
"""
    return ""


def snap_controls_str(snap_controls: dict) -> str:
    """
    Generate the snapControls section string for snappyHexMeshDict from a dictionary.

    Args:
        snap_controls (dict): Dictionary with the snapControls section settings.

    Returns:
        str: The snapControls section string.
    """
    return f"""
snapControls
{{
    nSmoothPatch {snap_controls.get('NSmoothPatch', 3)};
    tolerance {snap_controls.get('Tolerance', 2.0)};
    nSolveIter {snap_controls.get('NSolveIter', 30)};
    nRelaxIter {snap_controls.get('NRelaxIter', 5)};
    nFeatureSnapIter {snap_controls.get('NFeatureSnapIter', 10)};
}}
"""


def add_layers_controls_str(add_layers_controls: dict) -> str:
    """
    Generate the addLayersControls section string for snappyHexMeshDict from a
    dictionary.

    Args:
        add_layers_controls (dict): Dictionary with the addLayersControls section
            settings.

    Returns:
        str: The addLayersControls section string.
    """
    n_surface_layers = add_layers_controls.get('NSurfaceLayers', 0)
    return f"""
addLayersControls
{{
    relativeSizes {str(add_layers_controls.get('RelativeSizes', True)).lower()};
    layers
    {{
        airfoil
        {{
            nSurfaceLayers {n_surface_layers};
        }}
    }}
    expansionRatio {add_layers_controls.get('ExpansionRatio', 1.0)};
    finalLayerThickness {add_layers_controls.get('FinalLayerThickness', 0.3)};
    minThickness {add_layers_controls.get('MinThickness', 0.1)};
    nGrow {add_layers_controls.get('NGrow', 0)};
    featureAngle {add_layers_controls.get('FeatureAngle', 60)};
    nRelaxIter {add_layers_controls.get('NRelaxIter', 3)};
    nSmoothSurfaceNormals {add_layers_controls.get('NSmoothSurfaceNormals', 1)};
    nSmoothNormals {add_layers_controls.get('NSmoothNormals', 3)};
    nSmoothThickness {add_layers_controls.get('NSmoothThickness', 10)};
    maxFaceThicknessRatio {add_layers_controls.get('MaxFaceThicknessRatio', 0.5)};
    maxThicknessToMedialRatio {add_layers_controls.get('MaxThicknessToMedialRatio', 0.3)};
    minMedianAxisAngle {add_layers_controls.get('MinMedianAxisAngle', 90)};
    nBufferCellsNoExtrude {add_layers_controls.get('NBufferCellsNoExtrude', 0)};
    nLayerIter {add_layers_controls.get('NLayerIter', 50)};
}}
"""


def castellated_mesh_controls_str(
    chord: float,
    castellated_controls: dict,
    features_str: str,
    refinement_regions_str: str,
) -> str:
    """
    Generate the castellatedMeshControls section string for snappyHexMeshDict including
    refinementSurfaces creation with levels from the dict.

    Args:
        chord (float): Chord length of the airfoil.
        castellated_controls (dict): Dictionary with the castellatedMeshControls
            section settings.
        features_str (str): The features section string.
        refinement_regions_str (str): The refinementRegions section string.

    Returns:
        str: The castellatedMeshControls section string.
    """
    min_level = castellated_controls.get("MinSurfaceRefinementLevel", 0)
    max_level = castellated_controls.get("MaxSurfaceRefinementLevel", 0)
    refinement_surfaces_str = f"""    refinementSurfaces
    {{
        airfoil
        {{
            level ({min_level} {max_level});
        }}
    }}
"""

    return f"""
castellatedMeshControls
{{
    maxLocalCells {castellated_controls.get('MaxLocalCells', 100000)};
    maxGlobalCells {castellated_controls.get('MaxGlobalCells', 2000000)};
    minRefinementCells {castellated_controls.get('MinRefinementCells', 10)};
    maxLoadUnbalance {castellated_controls.get('MaxLoadUnbalance', 0.10)};
    nCellsBetweenLevels {castellated_controls.get('nCellsBetweenLevels', 3)};
{features_str}
{refinement_surfaces_str}
{refinement_regions_str}
    resolveFeatureAngle {castellated_controls.get('ResolveFeatureAngle', 30)};
    locationInMesh ({chord * 1.5} 0.0 0.0);
    allowFreeStandingZoneFaces true;
}}
"""


def mesh_quality_controls_str(mesh_quality_controls: dict) -> str:
    """
    Generate the meshQualityControls section string for snappyHexMeshDict from a
    dictionary.

    Args:
        mesh_quality_controls (dict): Dictionary with the meshQualityControls
            section settings.

    Returns:
        str: The meshQualityControls section string.
    """
    return f"""
meshQualityControls
{{
    maxNonOrtho {mesh_quality_controls.get('MaxNonOrtho', 65)};
    maxBoundarySkewness {mesh_quality_controls.get('MaxBoundarySkewness', 20)};
    maxInternalSkewness {mesh_quality_controls.get('MaxInternalSkewness', 4)};
    maxConcave {mesh_quality_controls.get('MaxConcave', 80)};
    minVol {mesh_quality_controls.get('MinVol', 1e-13)};
    minTetQuality {mesh_quality_controls.get('MinTetQuality', 1e-9)};
    minArea {mesh_quality_controls.get('MinArea', -1)};
    minTwist {mesh_quality_controls.get('MinTwist', 0.02)};
    minDeterminant {mesh_quality_controls.get('MinDeterminant', 0.001)};
    minFaceWeight {mesh_quality_controls.get('MinFaceWeight', 0.02)};
    minVolRatio {mesh_quality_controls.get('MinVolRatio', 0.01)};
    minTriangleTwist {mesh_quality_controls.get('MinTriangleTwist', -1)};
    nSmoothScale {mesh_quality_controls.get('NSmoothScale', 4)};
    errorReduction {mesh_quality_controls.get('ErrorReduction', 0.75)};
}}"""


def generate_cf_meshDict(
    thickness: float,
    location_in_mesh: tuple[float, float],
    max_cell_size: float,
    airfoil_ref_level: int,
    airfoil_ref_distance: float = 0.05
) -> str:
    """
    Generates a meshDict for cartesian2DMesh with airfoil surface refinement.

    Args:
        thickness: Mesh thickness in z-direction
        location_in_mesh: Point inside the fluid domain (x, y)
        max_cell_size: Maximum cell size in the domain
        airfoil_ref_level: Refinement level at the airfoil surface
        airfoil_ref_distance: Distance from airfoil surface where refinement is applied
    """
    refined_cell_size = max_cell_size / (2 ** airfoil_ref_level)

    return f"""FoamFile
{{
    version   2.0;
    format    ascii;
    class     dictionary;
    location  "system";
    object    meshDict;
}}

surfaceFile         "constant/triSurface/domain.fms";
maxCellSize         {max_cell_size};
locationInMesh      ({location_in_mesh[0]} {location_in_mesh[1]});

// Object-based refinements around the airfoil
localRefinement
{{
    airfoil
    {{
        additionalRefinementLevels 4;
    }}
}}
"""


def cf_mesh_dicts(airfoil: Airfoil, setup: Settings, system_dir: Path) -> None:
    """
    Generates a .fms file and a meshDict for the cartesian2DMesh workflow.
    Also saves individual STL files and a combined domain.stl for visualization.
    """
    case_dir = system_dir.parent
    tri_surface_dir = case_dir / "constant" / "triSurface"
    tri_surface_dir.mkdir(exist_ok=True, parents=True)

    cf = setup.mesh_settings.get("cfMesh", {})

    # --- 1. Define domain and patch names ---
    x_min, x_max, y_min, y_max = get_bounding_box(
        airfoil, setup.mesh_settings.get("BoundingBox", {}))
    thickness = cf.get("ZMax", 0.005) - cf.get("ZMin", -0.005)
    z_min = cf.get("ZMin", -0.005)
    z_max = cf.get("ZMax", 0.005)

    bc = setup.simulation_settings.get("BoundaryConditions", {})

    # This dictionary is for the edge-naming part of the .fms file
    edge_patch_names = {
        "inlet": bc.get("InletPatchName", "inlet"),
        "outlet": bc.get("OutletPatchName", "outlet"),
        "lowerWall": bc.get("LowerWallPatchName", "lowerWall"),
        "upperWall": bc.get("UpperWallPatchName", "upperWall"),
    }

    # --- 2. Create individual STL files for bounding box ---
    create_stl_bounding_box(
        x_min, x_max, y_min, y_max, z_min, z_max,
        edge_patch_names, tri_surface_dir
    )

    # --- 3. Create combined domain.stl by concatenating all STL files ---
    # The airfoil.stl should already exist from the main preparation script
    airfoil_stl_path = tri_surface_dir / "airfoil.stl"
    if not airfoil_stl_path.exists():
        raise FileNotFoundError(f"Airfoil STL not found at {airfoil_stl_path}")

    # List of all STL files to combine
    stl_files_to_combine = [airfoil_stl_path] + \
        [tri_surface_dir / f"{name}.stl" for name in edge_patch_names.values()]

    # Concatenate all STL files into domain.stl
    with open(tri_surface_dir / "domain.stl", "wb") as combined_file:
        for stl_file in stl_files_to_combine:
            if stl_file.exists():
                with open(stl_file, "rb") as individual_file:
                    combined_file.write(individual_file.read())

    # --- 4. Create the 2D domain .fms file ---
    x_contour = np.concatenate(
        (airfoil.upper_surface[0], airfoil.lower_surface[0][::-1]))
    y_contour = np.concatenate(
        (airfoil.upper_surface[1], airfoil.lower_surface[1][::-1]))

    export_domain_to_fms(
        x_contour, y_contour,
        x_min, x_max, y_min, y_max,
        patch_names=edge_patch_names,
        output_path=tri_surface_dir / "domain.fms"
    )

    # --- 5. Generate meshDict for cartesian2DMesh ---
    location_in_mesh = (x_min + 0.1 * airfoil.chord, 0.0)

    mesh_dict_content = generate_cf_meshDict(
        thickness=thickness,
        location_in_mesh=location_in_mesh,
        max_cell_size=cf.get("MaxCellSize", 0.05),
        airfoil_ref_level=cf.get('AirfoilRefinementMinLevel', 2)
    )
    with open(system_dir / "meshDict", "w") as f:
        f.write(mesh_dict_content)
