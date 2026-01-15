import numpy as np
from pathlib import Path
from templates.airfoil_template import Airfoil
from templates.initial_settings_template import Settings
from simulation.openfoam.block_mesh import get_bounding_box
from templates.openfoam_template_files.cf_mesh_files import generate_cf_mesh_dict
from utils.geometry import export_domain_to_fms, create_stl_bounding_box


def create_stl_domain(
    tri_surface_dir: Path,
    edge_patch_names: dict[str, str]
) -> None:
    """
    Create combined domain STL file from individual patch STL files.

    Args:
        tri_surface_dir: Directory containing STL files.
        edge_patch_names: Dictionary mapping edge types to patch names.

    Raises:
        FileNotFoundError: If airfoil.stl is not found.
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
    Generate .fms file and meshDict for cartesian2DMesh workflow.

    Creates individual STL files for boundary patches, a combined domain.stl
    file, and the meshDict configuration file.

    Args:
        airfoil: Airfoil object containing geometry data.
        setup: Settings object containing mesh parameters.
        system_dir: Path to the OpenFOAM case's system directory.
        stl_dir: Path where airfoil STL is located and where to save outputs.
                 Defaults to case_dir/constant/triSurface.
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
    mesh_dict_content = generate_cf_mesh_dict(
        airfoil=airfoil,
        location_in_mesh=location_in_mesh,
        max_cell_size=cf_mesh_settings.get("MaxCellSize", 0.02),
        local_refinement_dict=cf_mesh_settings.get("LocalRefinement", {}),
        boundary_layers_dict=cf_mesh_settings.get("BoundaryLayers", {}),
        object_refinement_dict=cf_mesh_settings.get("ObjectRefinements", {})
    )

    with open(system_dir / "meshDict", "w") as f:
        f.write(mesh_dict_content)
