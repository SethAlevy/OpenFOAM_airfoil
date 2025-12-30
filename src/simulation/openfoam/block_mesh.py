"""
BlockMesh generation utilities for OpenFOAM.
"""
from pathlib import Path
from airfoil.airfoil import Airfoil
from templates.initial_settings_template import Settings


def get_bounding_box(airfoil: Airfoil, bounding_box_settings: dict) -> tuple:
    """
    Calculate the bounding box dimensions for the mesh domain.

    Args:
        airfoil (Airfoil): Airfoil object containing geometry.
        bounding_box_settings (dict): Dictionary with bounding box parameters.

    Returns:
        tuple: Containing (x_min, x_max, y_min, y_max, z_min, z_max).
    """
    chord = airfoil.chord
    x_upstream = bounding_box_settings.get("InletDistance", 5.0) * chord
    x_downstream = bounding_box_settings.get("OutletDistance", 10.0) * chord
    y_top = bounding_box_settings.get("TopDistance", 5.0) * chord
    y_bottom = bounding_box_settings.get("BottomDistance", 5.0) * chord
    z_min = bounding_box_settings.get("ZMin", -0.5)
    z_max = bounding_box_settings.get("ZMax", 0.5)

    x_min = -x_upstream
    x_max = x_downstream
    y_min = -y_bottom
    y_max = y_top

    return x_min, x_max, y_min, y_max, z_min, z_max


def block_mesh_dict(
        airfoil: Airfoil,
        setup: Settings,
        output_path: Path
) -> None:
    """
    Generate blockMeshDict for OpenFOAM.

    Args:
        airfoil (Airfoil): Airfoil object containing geometry.
        setup (Settings): Settings object containing mesh parameters.
        output_path (Path): Path to save the blockMeshDict file.

    """
    mesh_settings = setup.mesh_settings
    block_mesh = mesh_settings.get("BlockMesh", {})

    x_min, x_max, y_min, y_max, z_min, z_max = get_bounding_box(
        airfoil, mesh_settings.get("BoundingBox", {}))

    nx = block_mesh.get("NX", 100)
    ny = block_mesh.get("NY", 50)
    nz = block_mesh.get("NZ", 1)

    scale = block_mesh.get("Scale", 1.0)

    bc = setup.simulation_settings.get("BoundaryConditions", {})
    inlet_name = bc.get("InletPatchName", "inlet")
    outlet_name = bc.get("OutletPatchName", "outlet")
    lower_wall_name = bc.get("LowerWallPatchName", "lowerWall")
    upper_wall_name = bc.get("UpperWallPatchName", "upperWall")

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale {scale};

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
    {inlet_name}
    {{
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }}
    {outlet_name}
    {{
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }}
    {lower_wall_name}
    {{
        type wall;
        faces
        (
            (0 1 5 4)
        );
    }}
    {upper_wall_name}
    {{
        type wall;
        faces
        (
            (3 7 6 2)
        );
    }}
    front
    {{
        type patch;
        faces
        (
            (4 5 6 7)
        );
    }}
    back
    {{
        type patch;
        faces
        (
            (0 3 2 1)
        );
    }}
);

mergePatchPairs
(
);
"""

    with open(output_path / "blockMeshDict", "w") as f:
        f.write(content)
