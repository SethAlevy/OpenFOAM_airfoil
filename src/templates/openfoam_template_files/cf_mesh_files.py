"""
cfMesh (cartesian2DMesh) dictionary template generation.
"""


def local_refinement(refinement_dict: dict) -> str:
    """
    Generate local refinement section for meshDict (patch-based).

    Args:
        refinement_dict: Dictionary with patch names as keys and refinement
            settings as values. Each should have 'CellSize'/'Level' and 'Thickness'.

    Returns:
        str: The localRefinement section, or empty string if no refinements.
    """
    if not refinement_dict:
        return ""

    refinement_entries = []
    for patch_name, settings in refinement_dict.items():
        thickness = settings.get("Thickness", 0.05)
        if "CellSize" in settings:
            cell_size = settings["CellSize"]
            entry = f"""    {patch_name}
    {{
        cellSize {cell_size};
        refinementThickness {thickness};
    }}"""
        else:
            level = settings.get("Level", 0)
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
        layers_dict: Dictionary with patch names as keys and layer settings
            as values. Each value should have: NLayers, ThicknessRatio,
            MaxFirstLayerThickness, AllowedDiscontinuity.

    Returns:
        str: The boundaryLayers section, or empty string if no layers.
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

    return f"""
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
    cell_size: float,
    name: str
) -> str:
    """
    Generate sphere refinement entry.

    Args:
        center: Center coordinates (x, y, z).
        radius: Sphere radius.
        cell_size: Desired cell size within sphere.
        name: Name identifier for this refinement object.

    Returns:
        str: Single sphere refinement entry.
    """
    return f"""    {name}
    {{
        type            sphere;
        centre          ({center[0]} {center[1]} {center[2]});
        radius          {radius};
        cellSize        {cell_size};
    }}
"""


def box_refinement(
    min_corner: tuple[float, float, float],
    max_corner: tuple[float, float, float],
    cell_size: float,
    name: str
) -> str:
    """
    Generate box refinement entry.

    Args:
        min_corner: Minimum corner (x_min, y_min, z_min).
        max_corner: Maximum corner (x_max, y_max, z_max).
        cell_size: Desired cell size within box.
        name: Name identifier for this refinement object.

    Returns:
        str: Single box refinement entry.
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
        centre          ({box_center[0]} {box_center[1]} {box_center[2]});
        lengthX         {x_length};
        lengthY         {y_length};
        lengthZ         1;
        cellSize        {cell_size};
    }}
"""


def object_refinements(airfoil, refinement_dict: dict) -> str:
    """
    Generate object refinements section for meshDict.

    Args:
        airfoil: Airfoil object containing geometry properties.
        refinement_dict: Dictionary with object names as keys and refinement
            settings as values. Each should specify Type ('sphere' or 'box')
            and type-specific parameters.

    Returns:
        str: The objectRefinements section, or empty string if no objects.
    """
    if not refinement_dict:
        return ""

    objects_dicts = ""
    for obj_name, obj_settings in refinement_dict.items():
        obj_type = obj_settings.get("Type", "").lower()

        if obj_type == "sphere":
            # Auto-detect special cases
            if "leadingedge" in obj_name.lower():
                centre = (0.0, 0.0, 0.0)
            elif "tip" in obj_name.lower():
                centre = (
                    airfoil.upper_surface[0][-1],
                    airfoil.upper_surface[1][-1],
                    0.0
                )
            else:
                centre = tuple(obj_settings.get("Center", (0.0, 0.0, 0.0)))

            objects_dicts += sphere_refinement(
                center=centre,
                radius=obj_settings.get("Radius", 0.2),
                cell_size=obj_settings.get("CellSize", 0.01),
                name=obj_name
            )

        elif obj_type == "box":
            x_min = obj_settings.get("XMin", -0.5)
            x_max = obj_settings.get("XMax", airfoil.chord * 2)
            y_min = obj_settings.get("YMin", -0.5)
            y_max = obj_settings.get("YMax", 0.5)

            objects_dicts += box_refinement(
                min_corner=tuple(obj_settings.get(
                    "MinCorner",
                    (x_min, y_min, 0.0)
                )),
                max_corner=tuple(obj_settings.get(
                    "MaxCorner",
                    (x_max, y_max, 0.0)
                )),
                cell_size=obj_settings.get("CellSize", 0.01),
                name=obj_name
            )

    return f"""
objectRefinements
{{
{objects_dicts}
}}
""" if objects_dicts else ""


def generate_cf_mesh_dict(
    airfoil,
    location_in_mesh: tuple[float, float],
    max_cell_size: float,
    local_refinement_dict: dict = None,
    boundary_layers_dict: dict = None,
    object_refinement_dict: dict = None
) -> str:
    """
    Generate complete meshDict for cartesian2DMesh.

    Args:
        airfoil: Airfoil object containing geometry data.
        location_in_mesh: Point inside the fluid domain (x, y).
        max_cell_size: Maximum cell size in the domain.
        local_refinement_dict: Patch-based refinement settings.
        boundary_layers_dict: Boundary layer settings per patch.
        object_refinement_dict: Object/region refinement settings.

    Returns:
        str: Complete meshDict file content.
    """
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
{object_section}{refinement_section}{layers_section}
"""
