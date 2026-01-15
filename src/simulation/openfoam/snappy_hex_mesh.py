from pathlib import Path
from templates.airfoil_template import Airfoil
from templates.initial_settings_template import Settings
from simulation.openfoam.block_mesh import get_bounding_box
from templates.openfoam_template_files.snappy_mex_mesh_file import (
    generate_snappy_hex_mesh_dict
)


def snappy_hex_mesh_dict(
    airfoil: Airfoil,
    setup: Settings,
    output_path: Path
) -> None:
    """
    Generate snappyHexMeshDict with refinement, layers, and defaults aligned to cfMesh
    utilities.

    Args:
        airfoil: Airfoil object with geometry data
        setup: Settings object containing mesh and simulation configuration
        output_path: Path to the output directory
    """
    mesh_settings = setup.mesh_settings
    snappy = mesh_settings.get("SnappyHexMesh", {})

    x_min, x_max, y_min, y_max, z_min, z_max = get_bounding_box(
        airfoil, mesh_settings.get("BoundingBox", {})
    )

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
