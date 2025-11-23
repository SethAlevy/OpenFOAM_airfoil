from templates.airfoil_template import Airfoil
from templates.initial_settings_template import Settings
from pathlib import Path


def block_mesh_dict(airfoil: Airfoil, setup: Settings, output_path: Path) -> None:
    """
    Fill the blockMeshDict templates according to the selected simulation setup and
    save the file in the given path.

    Args:
        airfoil (Airfoil): Object with the generated airfoil geometry.
        setup (Settings): Object with the simulation settings.
        output_path (Path): Path where the blockMeshDict file will be saved.
    """
    bounding_box = setup.mesh_settings.get("BoundingBox")
    x_min, x_max, y_min, y_max = get_bounding_box(airfoil, bounding_box)
    cell_size = bounding_box["BaseCellSize"]

    content = generate_block_mesh_dict(
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        nx=int((x_max - x_min) / cell_size),
        ny=int((y_max - y_min) / cell_size),
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
    z_min: float = 0,
    z_max: float = 0.01,
    nz: int = 1
):
    """
    Generate a simple BlockMesh dict for a 2D airfoil case.

    Args:
        x_min (float): Minimum x-coordinate of the domain.
        x_max (float): Maximum x-coordinate of the domain.
        y_min (float): Minimum y-coordinate of the domain.
        y_max (float): Maximum y-coordinate of the domain.
        nx (int): Number of cells in the x-direction.
        ny (int): Number of cells in the y-direction.
        z_min (float, optional): Minimum z-coordinate of the domain. Defaults to 0.
        z_max (float, optional): Maximum z-coordinate of the domain. Defaults to 0.01.
        nz (int, optional): Number of cells in the z-direction. Defaults to 1
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
    inlet
    {{
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }}
    outlet
    {{
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }}
    lowerWall
    {{
        type wall;
        faces
        (
            (0 1 5 4)
        );
    }}
    upperWall
    {{
        type wall;
        faces
        (
            (3 7 6 2)
        );
    }}
    frontAndBack
    {{
        type empty;
        faces
        (
            (0 3 2 1)
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
        setup (Settings): Object with the simulation settings.
        output_path (Path): Path where the snappyHexMeshDict file will be saved.
    """
    snappy_hex_mesh_settings = setup.mesh_settings.get("SnappyHexMesh")
    add_layers = snappy_hex_mesh_settings.get("AddLayers", False)
    n_layers = snappy_hex_mesh_settings.get("NSurfaceLayers", 0)
    min_surface_refinement = snappy_hex_mesh_settings.get("MinSurfaceRefinementLevel", 0)
    max_surface_refinement = snappy_hex_mesh_settings.get("MaxSurfaceRefinementLevel", 0)
    feature_refinement_level = snappy_hex_mesh_settings.get("FeatureRefinementLevel", 0)
    refinement_box_level = snappy_hex_mesh_settings.get("RefinementBoxLevel", 0)
    sphere_tip_level = snappy_hex_mesh_settings.get("SphereTipLevel", 0)
    sphere_tip_radius_scale = snappy_hex_mesh_settings.get("SphereTipRadius", 0)
    sphere_leading_level = snappy_hex_mesh_settings.get("SphereLeadingEdgeLevel", 0)
    sphere_leading_radius_scale = snappy_hex_mesh_settings.get("SphereLeadingEdgeRadius", 0)

    if refinement_box_level > 0:
        refinement_box = refinement_box_dict(airfoil)
    else:
        refinement_box = ''
    if sphere_tip_level > 0:
        sphere_tip = refine_sphere_dict(airfoil, "tip", sphere_tip_radius_scale)
    else:
        sphere_tip = ''
    if sphere_leading_level > 0:
        sphere_leading = refine_sphere_dict(airfoil, "leading", sphere_leading_radius_scale)
    else:
        sphere_leading = ''
    content = generate_snappy_hex_mesh_dict(
        add_layers=add_layers,
        refinement_box=refinement_box,
        refinement_box_level=refinement_box_level,
        sphere_tip=sphere_tip,
        sphere_tip_level=sphere_tip_level,
        sphere_leading=sphere_leading,
        sphere_leading_level=sphere_leading_level,
        surface_refinement_levels=(min_surface_refinement, max_surface_refinement),
        n_surface_layers=n_layers,
        feature_refinement_level=feature_refinement_level,
    )

    output_path.mkdir(parents=True, exist_ok=True)
    with open(output_path / "snappyHexMeshDict", 'w') as file:
        file.write(content)


def refinement_box_dict(Airfoil: Airfoil) -> dict:
    chord = Airfoil.chord
    upper_surface = Airfoil.upper_surface
    lower_surface = Airfoil.lower_surface

    return {
        "min": [min(lower_surface[0]) - chord, min(lower_surface[1]) - chord, -0.1],
        "max": [max(upper_surface[0]) + chord, max(upper_surface[1]) + chord, 0.1],
    }


def refine_sphere_dict(Airfoil: Airfoil, position: str, radius_scale: float) -> dict:
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


def generate_snappy_hex_mesh_dict(
    add_layers: bool = True,
    refinement_box: dict = '',
    refinement_box_level: int = 0,
    sphere_tip: dict = '',
    sphere_tip_level: int = 0,
    sphere_leading: dict = '',
    sphere_leading_level: int = 0,
    surface_refinement_levels: tuple = (0, 0),
    n_surface_layers: int = 0,
    feature_refinement_level: int = 0,
):
    """
    Generate a basic snappyHexMeshDict for a 2D airfoil case, with optional refinement box and spheres.

    Args:
        add_layers (bool): Enable layer addition.
        refinement_box (dict, optional): Dict with 'min' and 'max' keys for refinement box.
        sphere_tip (dict, optional): Dict with 'center' (list of 3 floats) and 'radius' (float).
        sphere_leading (dict, optional): Dict with 'center' (list of 3 floats) and 'radius' (float).
        refinement_levels (tuple): Refinement levels for airfoil surface (e.g., (3, 4)).
        n_surface_layers (int): Number of surface layers for airfoil.
        feature_level (int, optional): Refinement level for the feature.
        output_path (str): Path to save the snappyHexMeshDict.
    """
    ref_box_geo, ref_box_reg = box_geometry_and_refinement_str(
        "refinementBox", refinement_box, refinement_box_level)
    tip_sphere_geo, tip_sphere_reg = sphere_geometry_and_refinement_str(
        "tipSphere", sphere_tip, sphere_tip_level)
    lead_sphere_geo, lead_sphere_reg = sphere_geometry_and_refinement_str(
        "leadingEdgeSphere", sphere_leading, sphere_leading_level)

    # Features section
    features_str = feature_str(feature_refinement_level)

    # Geometry section
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

    # Refinement regions section
    refinement_regions_str = "    refinementRegions\n    {\n"
    refinement_regions_str += ref_box_reg
    refinement_regions_str += tip_sphere_reg
    refinement_regions_str += lead_sphere_reg
    refinement_regions_str += "    }\n"

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}

castellatedMesh true;
snap            true;
addLayers       {str(add_layers).lower()};

geometry
{{
{geometry_section}
}}

castellatedMeshControls
{{
    maxLocalCells 100000;
    maxGlobalCells 2000000;
    minRefinementCells 10;
    nCellsBetweenLevels 3;

{features_str}    refinementSurfaces
    {{
        airfoil
        {{
            level ({surface_refinement_levels[0]} {surface_refinement_levels[1]});
        }}
    }}

{refinement_regions_str}
    resolveFeatureAngle 30;
    locationInMesh (0 0 0.005);
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 30;
    nRelaxIter 5;
}}

addLayersControls
{{
    relativeSizes true;
    layers
    {{
        airfoil
        {{
            nSurfaceLayers {n_surface_layers};
        }}
    }}
    expansionRatio 1.0;
    finalLayerThickness 0.3;
    minThickness 0.1;
    nGrow 0;
    featureAngle 60;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}}

meshQualityControls
{{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-13;
    minTetQuality 1e-9;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.02;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
}}

debug 0;
mergeTolerance 1E-6;
"""
    return content


def box_geometry_and_refinement_str(name, box, level):
    """
    Generate geometry and refinement region strings for a box.
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


def sphere_geometry_and_refinement_str(name, sphere, level):
    """
    Generate geometry and refinement region strings for a sphere.
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
    Generate the features section string for snappyHexMeshDict.
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
