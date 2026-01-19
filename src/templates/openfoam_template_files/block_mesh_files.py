"""
blockMesh dictionary template generation.

This module contains functions to generate OpenFOAM blockMeshDict files
for structured hex meshing.
"""


def vertices_section(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    z_min: float,
    z_max: float
) -> str:
    """
    Generate the vertices section of blockMeshDict.

    Creates 8 vertices for a rectangular hex block:
    0: (x_min, y_min, z_min)
    1: (x_max, y_min, z_min)
    2: (x_max, y_max, z_min)
    3: (x_min, y_max, z_min)
    4: (x_min, y_min, z_max)
    5: (x_max, y_min, z_max)
    6: (x_max, y_max, z_max)
    7: (x_min, y_max, z_max)

    Args:
        x_min (float): Domain lower bound in x direction.
        x_max (float): Domain upper bound in x direction.
        y_min (float): Domain lower bound in y direction.
        y_max (float): Domain upper bound in y direction.
        z_min (float): Domain lower bound in z direction.
        z_max (float): Domain upper bound in z direction.

    Returns:
        str: Formatted vertices section.
    """
    return f"""vertices
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
"""


def blocks_section(
    nx: int,
    ny: int,
    nz: int
) -> str:
    """
    Generate the blocks section of blockMeshDict.

    Creates a single hex block with specified cell divisions.

    Args:
        nx (int): Number of cells in x direction.
        ny (int): Number of cells in y direction.
        nz (int): Number of cells in z direction.

    Returns:
        str: Formatted blocks section.
    """
    return f"""blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
);
"""


def boundary_section(
    inlet_patch: str,
    outlet_patch: str,
    lower_wall_patch: str,
    upper_wall_patch: str,
    front_patch: str = "front",
    back_patch: str = "back"
) -> str:
    """
    Generate the boundary section of blockMeshDict.

    Creates 6 boundary patches (inlet, outlet, lower wall, upper wall,
    front/back symmetry planes).

    Args:
        inlet_patch (str): Name of inlet boundary patch.
        outlet_patch (str): Name of outlet boundary patch.
        lower_wall_patch (str): Name of lower wall patch.
        upper_wall_patch (str): Name of upper wall patch.
        front_patch (str): Name of front symmetry plane.
        back_patch (str): Name of back symmetry plane.

    Returns:
        str: Formatted boundary section.
    """
    return f"""boundary
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
        type wall;
        faces
        (
            (0 1 5 4)
        );
    }}
    {upper_wall_patch}
    {{
        type wall;
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
            (4 5 6 7)
        );
    }}
    {back_patch}
    {{
        type symmetryPlane;
        faces
        (
            (0 3 2 1)
        );
    }}
);
"""


def generate_block_mesh_dict(
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    z_min: float,
    z_max: float,
    nx: int,
    ny: int,
    nz: int,
    scale: float = 1.0,
    inlet_patch: str = "inlet",
    outlet_patch: str = "outlet",
    lower_wall_patch: str = "lowerWall",
    upper_wall_patch: str = "upperWall",
    front_patch: str = "front",
    back_patch: str = "back"
) -> str:
    """
    Generate complete blockMeshDict file content.

    Args:
        x_min (float): Domain lower bound in x direction.
        x_max (float): Domain upper bound in x direction.
        y_min (float): Domain lower bound in y direction.
        y_max (float): Domain upper bound in y direction.
        z_min (float): Domain lower bound in z direction.
        z_max (float): Domain upper bound in z direction.
        nx (int): Number of cells in x direction.
        ny (int): Number of cells in y direction.
        nz (int): Number of cells in z direction.
        scale (float): Scale factor for all coordinates.
        inlet_patch (str): Name of inlet boundary patch.
        outlet_patch (str): Name of outlet boundary patch.
        lower_wall_patch (str): Name of lower wall patch.
        upper_wall_patch (str): Name of upper wall patch.
        front_patch (str): Name of front symmetry plane.
        back_patch (str): Name of back symmetry plane.

    Returns:
        str: Complete blockMeshDict file content.
    """
    vertices = vertices_section(x_min, x_max, y_min, y_max, z_min, z_max)
    blocks = blocks_section(nx, ny, nz)
    boundary = boundary_section(
        inlet_patch, outlet_patch, lower_wall_patch, upper_wall_patch,
        front_patch, back_patch
    )

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

scale {scale};

{vertices}
{blocks}
edges
(
);

{boundary}
mergePatchPairs
(
);
"""

    return content
