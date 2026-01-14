from pathlib import Path
from templates.airfoil_template import Airfoil
from templates.initial_settings_template import Settings
from simulation.openfoam.block_mesh import get_bounding_box


def _bool_to_foam(value: bool) -> str:
    """Return OpenFOAM-friendly boolean text."""
    return "true" if value else "false"


def _level_pair(raw_level, default: tuple[int, int]) -> tuple[int, int]:
    """Normalize refinement level input to a pair."""
    if raw_level is None:
        return default
    if isinstance(raw_level, (int, float)):
        return (int(raw_level), int(raw_level))
    if isinstance(raw_level, (list, tuple)):
        if len(raw_level) == 1:
            return (int(raw_level[0]), int(raw_level[0]))
        return (int(raw_level[0]), int(raw_level[1]))
    return default


def _refinement_surfaces(
    surface_dict: dict,
    default_patch: str,
    default_levels: tuple[int, int]
) -> str:
    """Build the refinementSurfaces section."""
    if not surface_dict:
        surface_dict = {default_patch: {"Level": default_levels}}

    surface_entries = []
    for patch_name, settings in surface_dict.items():
        level_min, level_max = _level_pair(settings.get("Level"), default_levels)
        patch_type = settings.get("PatchType")
        face_zone = settings.get("FaceZone")

        extra = []
        if patch_type:
            extra.append(
                """            patchInfo
            {{
                type {patch_type};
            }}""".format(patch_type=patch_type)
            )
        if face_zone:
            extra.append(f"            faceZone {face_zone};")

        extra_block = "\n" + "\n".join(extra) if extra else ""

        surface_entries.append(
            f"""        {patch_name}
        {{
            level ({level_min} {level_max});{extra_block}
        }}"""
        )

    joined = "\n\n".join(surface_entries)
    return f"""    refinementSurfaces
    {{
{joined}
    }}"""


def _box_geometry(min_corner: tuple[float, float, float], max_corner: tuple[float, float, float], name: str) -> str:
    """Create searchableBox geometry entry."""
    return f"""    {name}
    {{
        type searchableBox;
        min ({min_corner[0]} {min_corner[1]} {min_corner[2]});
        max ({max_corner[0]} {max_corner[1]} {max_corner[2]});
    }}"""


def _sphere_geometry(center: tuple[float, float, float], radius: float, name: str) -> str:
    """Create searchableSphere geometry entry."""
    return f"""    {name}
    {{
        type searchableSphere;
        centre ({center[0]} {center[1]} {center[2]});
        radius {radius};
    }}"""


def _refinement_regions(
    airfoil: Airfoil,
    region_dict: dict,
    default_levels: tuple[int, int]
) -> tuple[str, str]:
    """Build geometry snippets and refinementRegions section from object refinements."""
    if not region_dict:
        return "", "    refinementRegions\n    {\n    }"

    geometry_entries: list[str] = []
    region_entries: list[str] = []

    for name, settings in region_dict.items():
        ref_type = settings.get("Type", "box").lower()
        level_min, level_max = _level_pair(settings.get("Level"), default_levels)
        region_level = max(level_min, level_max)
        mode = settings.get("Mode", "inside")

        if ref_type == "sphere":
            if "leadingedge" in name.lower():
                center = (0.0, 0.0, 0.0)
            elif "tip" in name.lower():
                center = (
                    airfoil.upper_surface[0][-1],
                    airfoil.upper_surface[1][-1],
                    0.0,
                )
            else:
                center = tuple(settings.get("Center", (0.0, 0.0, 0.0)))
            radius = settings.get("Radius", 0.1)
            geometry_entries.append(_sphere_geometry(center, radius, name))
        else:
            x_min = settings.get("XMin", -0.5)
            x_max = settings.get("XMax", airfoil.chord * 2)
            y_min = settings.get("YMin", -0.5)
            y_max = settings.get("YMax", 0.5)
            z_min = settings.get("ZMin", -0.001 * airfoil.chord)
            z_max = settings.get("ZMax", 0.001 * airfoil.chord)
            min_corner = tuple(settings.get("MinCorner", (x_min, y_min, z_min)))
            max_corner = tuple(settings.get("MaxCorner", (x_max, y_max, z_max)))
            geometry_entries.append(_box_geometry(min_corner, max_corner, name))

        region_entries.append(
            f"""        {name}
        {{
            mode {mode};
            levels ((1E15 {region_level}));
        }}"""
        )

    geometry_block = "\n" + "\n\n".join(geometry_entries) + "\n"
    regions_block = "\n".join(region_entries)
    return geometry_block, f"""    refinementRegions
    {{
{regions_block}
    }}"""


def _feature_section(level: int, file_name: str = "airfoil.eMesh") -> str:
    """Build the features section if requested."""
    if level and level > 0:
        return f"""    features
    (
        {{
            file "{file_name}";
            level {level};
        }}
    );"""
    return "    features\n    (\n    );"


def _layers_section(layer_dict: dict, add_layers_settings: dict) -> tuple[str, bool]:
    """Build layers subsection and flag whether layers are enabled."""
    layer_entries = []
    for patch_name, settings in layer_dict.items():
        n_layers = settings.get("NLayers", 0)
        if n_layers and n_layers > 0:
            layer_entries.append(
                f"""        {patch_name}
        {{
            nSurfaceLayers {n_layers};
        }}"""
            )

    expansion_ratio = add_layers_settings.get("ExpansionRatio", 1.15)
    final_layer_thickness = add_layers_settings.get("FinalLayerThickness", 0.25)
    min_thickness = add_layers_settings.get("MinThickness", 0.05)
    n_grow = add_layers_settings.get("NGrow", 2)
    feature_angle = add_layers_settings.get("FeatureAngle", 150)
    n_relax_iter = add_layers_settings.get("NRelaxIter", 5)
    n_smooth_surface_normals = add_layers_settings.get("NSmoothSurfaceNormals", 3)
    n_smooth_normals = add_layers_settings.get("NSmoothNormals", 10)
    n_smooth_thickness = add_layers_settings.get("NSmoothThickness", 10)
    max_face_thickness_ratio = add_layers_settings.get("MaxFaceThicknessRatio", 0.6)
    max_thickness_to_medial_ratio = add_layers_settings.get(
        "MaxThicknessToMedialRatio", 0.9)
    min_medial_axis_angle = add_layers_settings.get("MinMedialAxisAngle", 50)
    n_buffer_cells_no_extrude = add_layers_settings.get("NBufferCellsNoExtrude", 2)
    n_layer_iter = add_layers_settings.get("NLayerIter", 150)

    layers_block = "\n".join(layer_entries)
    relative_sizes = add_layers_settings.get("RelativeSizes", True)
    add_layers_enabled = bool(layer_entries)

    section = f"""addLayersControls
{{
    relativeSizes {_bool_to_foam(relative_sizes)};
    layers
    {{
{layers_block}
    }}
    expansionRatio {expansion_ratio};
    finalLayerThickness {final_layer_thickness};
    minThickness {min_thickness};
    nGrow {n_grow};
    featureAngle {feature_angle};
    nRelaxIter {n_relax_iter};
    nSmoothSurfaceNormals {n_smooth_surface_normals};
    nSmoothNormals {n_smooth_normals};
    nSmoothThickness {n_smooth_thickness};
    maxFaceThicknessRatio {max_face_thickness_ratio};
    maxThicknessToMedialRatio {max_thickness_to_medial_ratio};
    minMedialAxisAngle {min_medial_axis_angle};
    nBufferCellsNoExtrude {n_buffer_cells_no_extrude};
    nLayerIter {n_layer_iter};
}}"""

    return section, add_layers_enabled


def generate_snappy_hex_mesh_dict(
    airfoil: Airfoil,
    setup: Settings,
    location_in_mesh: tuple[float, float, float],
    snappy_settings: dict,
    surface_refinement_dict: dict,
    region_refinement_dict: dict,
    boundary_layers_dict: dict
) -> str:
    """Compose the full snappyHexMeshDict content."""

    default_levels = (
        snappy_settings.get("RefinementMinLevel", 2),
        snappy_settings.get("RefinementMaxLevel", 4),
    )

    # Geometry and refinement regions
    geometry_extra, refinement_regions = _refinement_regions(
        airfoil=airfoil,
        region_dict=region_refinement_dict,
        default_levels=default_levels,
    )

    geometry_section = f"""geometry
{{
    airfoil.stl
    {{
        type triSurfaceMesh;
        name airfoil;
    }}{geometry_extra}
}};"""

    # Castellated mesh controls
    castellated_defaults = {
        "MaxLocalCells": 2000000,
        "MaxGlobalCells": 8000000,
        "MinRefinementCells": 0,
        "MaxLoadUnbalance": 0.10,
        "NCellsBetweenLevels": 3,
        "ResolveFeatureAngle": 30,
        "AllowFreeStandingZoneFaces": True,
        "FeatureRefinementLevel": 0,
    }

    castellated = {
        **castellated_defaults,
        **snappy_settings.get("CastellatedMeshControls", {}),
    }

    features_section = _feature_section(
        level=castellated.get("FeatureRefinementLevel", 0)
    )
    refinement_surfaces = _refinement_surfaces(
        surface_dict=surface_refinement_dict,
        default_patch=setup.simulation_settings
        .get("BoundaryConditions", {})
        .get("AirfoilPatchName", "airfoil"),
        default_levels=default_levels,
    )

    # Snap controls
    snap_defaults = {
        "NSmoothPatch": 3,
        "Tolerance": 2.0,
        "NSolveIter": 30,
        "NRelaxIter": 5,
    }
    snap_controls = {**snap_defaults, **snappy_settings.get("SnapControls", {})}

    # Layer controls
    add_layers_section, add_layers_enabled = _layers_section(
        boundary_layers_dict,
        snappy_settings.get("AddLayersControls", {}),
    )

    # Mesh quality controls
    quality_defaults = {
        "MaxNonOrtho": 65,
        "MaxBoundarySkewness": 20,
        "MaxInternalSkewness": 4,
        "MaxConcave": 80,
        "MinVol": 1e-13,
        "MinTetQuality": 1e-30,
        "MinArea": -1,
        "MinTwist": 0.02,
        "MinDeterminant": 0.001,
        "MinFaceWeight": 0.02,
        "MinVolRatio": 0.01,
        "MinTriangleTwist": -1,
        "NSmoothScale": 4,
        "ErrorReduction": 0.75,
    }
    relaxed_defaults = {
        "MaxNonOrtho": 80,
        "MaxBoundarySkewness": 25,
        "MaxInternalSkewness": 4.5,
        "MaxConcave": 85,
        "MinVol": 1e-15,
        "MinTetQuality": 1e-35,
        "MinTwist": 0.005,
        "MinDeterminant": 0.0005,
        "MinFaceWeight": 0.01,
        "MinVolRatio": 0.005,
    }
    mq_user = snappy_settings.get("MeshQualityControls", {})
    mesh_quality = {**quality_defaults, **{k: v for k, v in mq_user.items() if k not in ("relaxed", "Relaxed")}}
    relaxed = {**relaxed_defaults, **mq_user.get("relaxed", mq_user.get("Relaxed", {}))}

    debug = snappy_settings.get("Debug", 0)
    merge_tolerance = snappy_settings.get("MergeTolerance", 1e-6)

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}

castellatedMesh true;
snap            true;
addLayers       {_bool_to_foam(add_layers_enabled)};

{geometry_section}

castellatedMeshControls
{{
    maxLocalCells {castellated['MaxLocalCells']};
    maxGlobalCells {castellated['MaxGlobalCells']};
    minRefinementCells {castellated['MinRefinementCells']};
    maxLoadUnbalance {castellated['MaxLoadUnbalance']};
    nCellsBetweenLevels {castellated['NCellsBetweenLevels']};

{features_section}

{refinement_surfaces}

    resolveFeatureAngle {castellated['ResolveFeatureAngle']};

{refinement_regions}

    locationInMesh ({location_in_mesh[0]} {location_in_mesh[1]} {location_in_mesh[2]});
    allowFreeStandingZoneFaces {_bool_to_foam(castellated['AllowFreeStandingZoneFaces'])};
}}

snapControls
{{
    nSmoothPatch {snap_controls['NSmoothPatch']};
    tolerance {snap_controls['Tolerance']};
    nSolveIter {snap_controls['NSolveIter']};
    nRelaxIter {snap_controls['NRelaxIter']};
}}

{add_layers_section}

meshQualityControls
{{
    maxNonOrtho {mesh_quality['MaxNonOrtho']};
    maxBoundarySkewness {mesh_quality['MaxBoundarySkewness']};
    maxInternalSkewness {mesh_quality['MaxInternalSkewness']};
    maxConcave {mesh_quality['MaxConcave']};
    minVol {mesh_quality['MinVol']};
    minTetQuality {mesh_quality['MinTetQuality']};
    minArea {mesh_quality['MinArea']};
    minTwist {mesh_quality['MinTwist']};
    minDeterminant {mesh_quality['MinDeterminant']};
    minFaceWeight {mesh_quality['MinFaceWeight']};
    minVolRatio {mesh_quality['MinVolRatio']};
    minTriangleTwist {mesh_quality['MinTriangleTwist']};

    nSmoothScale {mesh_quality['NSmoothScale']};
    errorReduction {mesh_quality['ErrorReduction']};

    relaxed
    {{
        maxNonOrtho {relaxed['MaxNonOrtho']};
        maxBoundarySkewness {relaxed['MaxBoundarySkewness']};
        maxInternalSkewness {relaxed['MaxInternalSkewness']};
        maxConcave {relaxed['MaxConcave']};
        minVol {relaxed['MinVol']};
        minTetQuality {relaxed['MinTetQuality']};
        minTwist {relaxed['MinTwist']};
        minDeterminant {relaxed['MinDeterminant']};
        minFaceWeight {relaxed['MinFaceWeight']};
        minVolRatio {relaxed['MinVolRatio']};
    }}
}}

debug {debug};
mergeTolerance {merge_tolerance};
"""

    return content


def snappy_hex_mesh_dict(
    airfoil: Airfoil,
    setup: Settings,
    output_path: Path
) -> None:
    """Generate snappyHexMeshDict with refinement, layers, and defaults aligned to cfMesh utilities."""
    mesh_settings = setup.mesh_settings
    snappy = mesh_settings.get("SnappyHexMesh", {})

    x_min, x_max, y_min, y_max, z_min, z_max = get_bounding_box(
        airfoil, mesh_settings.get("BoundingBox", {})
    )

    # Pick a location safely inside the domain: offset from the STL plane
    default_location = (
        -0.5 * airfoil.chord,
        0.0,
        0.02,   # small z-offset
    )
    location_in_mesh = tuple(snappy.get("LocationInMesh", default_location))

    surface_refinement = snappy.get(
        "LocalRefinement", snappy.get("SurfaceRefinement", {}))
    region_refinement = snappy.get("ObjectRefinements", {})
    boundary_layers = snappy.get("BoundaryLayers", {})

    content = generate_snappy_hex_mesh_dict(
        airfoil=airfoil,
        setup=setup,
        location_in_mesh=location_in_mesh,
        snappy_settings=snappy,
        surface_refinement_dict=surface_refinement,
        region_refinement_dict=region_refinement,
        boundary_layers_dict=boundary_layers,
    )

    with open(output_path / "snappyHexMeshDict", "w") as file:
        file.write(content)
