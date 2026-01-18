from pathlib import Path
from simulation.openfoam.block_mesh import block_mesh_dict
from simulation.openfoam.snappy_hex_mesh import snappy_hex_mesh_dict
from simulation.openfoam.cf_mesh_2d import cf_mesh_dict
from simulation.openfoam.system_dir import (
    add_force_coeffs_dict_from_bc,
    control_dict,
    fv_schemes_dict,
    surface_feature_extract_dict,
    decompose_par_dict,
    fv_solution_dict,
    extrude_mesh_dict
)
from simulation.openfoam.constant_dir import transport_properties_dict, \
    turbulence_properties_dict
from simulation.openfoam.boundary_condition import BoundaryConditions
from templates.python_template_files.airfoil_template import Airfoil
from templates.python_template_files.initial_settings_template import Settings
from utils.logger import SimpleLogger


def prepare_openfoam_case(
        working_path: Path,
        case_name: str,
        airfoil: Airfoil,
        setup: Settings
) -> None:
    """
    Prepares the OpenFOAM case by setting up necessary files and directories.

    Args:
        working_path (Path): The base working directory path.
        case_name (str): The name of the case.
        airfoil (Airfoil): The airfoil object containing geometry data.
        setup (Settings): The simulation settings.
    """
    create_directory_structure(working_path, case_name)

    mesher = setup.mesh_settings.get("Mesher")
    z_min = setup.mesh_settings.get("BoundingBox", {}).get("ZMin", -0.5)
    z_max = setup.mesh_settings.get("BoundingBox", {}).get("ZMax", 0.5)

    if mesher and mesher.lower() in ["placeholder2d"]:
        airfoil.to_stl(
            output_path=(
                working_path / case_name / "constant" / "triSurface" / "airfoil.stl"
            ),
            dimension=2
        )
    elif mesher and mesher.lower() in ["cfmesh", "snappyhexmesh"]:
        thickness = (z_max - z_min) * 1.01
        airfoil.to_stl(
            output_path=(
                working_path / case_name / "constant" / "triSurface" / "airfoil.stl"
            ),
            thickness=thickness,
            dimension=3
        )

    create_meshing_files(working_path / case_name, airfoil, setup, mesher)
    create_system_files(working_path / case_name, setup)
    create_constant_files(working_path / case_name, setup)
    bc = create_boundary_conditions_files(working_path / case_name / "0", setup)
    add_force_coeffs_dict_from_bc(
        working_path / case_name / "system" / "controlDict",
        bc,
        chord=airfoil.chord,
        span=abs(z_max - z_min),
        angle_of_attack_deg=airfoil.alpha
    )


def create_directory_structure(working_path: Path, case_name: str) -> None:
    """
    Creates the required directory structure for the OpenFOAM case.

    Args:
        working_path (Path): The base working directory path.
        case_name (str): The name of the case.
    """
    working_path.mkdir(parents=True, exist_ok=True)
    (working_path / case_name).mkdir(exist_ok=True)

    (working_path / case_name / "0").mkdir(exist_ok=True)
    (working_path / case_name / "constant").mkdir(exist_ok=True)
    (working_path / case_name / "constant" / "triSurface").mkdir(exist_ok=True)
    (working_path / case_name / "system").mkdir(exist_ok=True)


def create_meshing_files(
        case_path: Path,
        airfoil: Airfoil,
        setup: Settings,
        mesher: str
) -> None:
    """
    Creates meshing files for the OpenFOAM case.

    Args:
        case_path (Path): The path to the case directory.
        airfoil (Airfoil): The airfoil object containing geometry data.
        setup (Settings): The simulation settings.
        mesher (str): The meshing method to use (e.g., "cfMesh", "snappyHexMesh").
    """

    block_mesh_dict(airfoil, setup, (case_path / "system"))

    if mesher and mesher.lower() == "snappyhexmesh":
        SimpleLogger.log("SnappyHexMesh selected as mesher. Creating file.")
        snappy_hex_mesh_dict(airfoil, setup, (case_path / "system"))
    elif mesher and mesher.lower() == "cfmesh":
        SimpleLogger.log("cfMesh selected as mesher. No additional files needed.")
        cf_mesh_dict(airfoil, setup, (case_path / "system"))
    else:
        SimpleLogger.log("No specific mesher selected or unrecognized mesher.")


def create_system_files(case_path: Path, setup: Settings) -> None:
    """
    Creates system files for the OpenFOAM case.

    Args:
        case_path (Path): The path to the case directory.
        setup (Settings): The simulation settings.
    """
    control_dict(setup, (case_path / "system" / "controlDict"))
    fv_solution_dict(setup, (case_path / "system" / "fvSolution"))
    fv_schemes_dict(setup, (case_path / "system" / "fvSchemes"))

    n_processors = (
        setup.simulation_settings
        .get("Decomposition", {})
        .get("NumberOfSubdomains", 0)
    )

    if n_processors and n_processors > 1:
        decompose_par_dict(setup, (case_path / "system" / "decomposeParDict"))

    mesher = setup.mesh_settings.get("Mesher", "").lower()

    if mesher == "snappyhexmesh":
        feature_refinement = (
            setup.mesh_settings
            .get("SnappyHexMesh", {})
            .get("CastellatedMeshControls", {})
            .get("FeatureRefinementLevel", 0)
        )
        if feature_refinement and feature_refinement > 0:
            surface_feature_extract_dict(
                setup,
                (case_path / "system" / "surfaceFeatureExtractDict")
            )

        extrude_mesh_dict(setup, (case_path / "system" / "extrudeMeshDict"))


def create_constant_files(case_path: Path, setup: Settings) -> None:
    """
    Creates constant directory files for the OpenFOAM case.

    Args:
        case_path (Path): The path to the case directory.
        setup (Settings): The simulation settings.
    """
    transport_properties_dict(setup, (case_path / "constant" / "transportProperties"))
    turbulence_properties_dict(setup, (case_path / "constant" / "turbulenceProperties"))


def create_boundary_conditions_files(bc_path: Path, setup: Settings) -> None:
    """
    Creates boundary condition files for the OpenFOAM case.

    Args:
        bc_path (Path): The path to the boundary conditions directory ("0").
        setup (Settings): The simulation settings.
    """
    bc = BoundaryConditions(setup=setup)
    velocity_content = bc.velocity_bc()
    pressure_content = bc.pressure_bc()

    bc.write_bc(velocity_content, bc_path / "U")
    bc.write_bc(pressure_content, bc_path / "p")

    turbulence_content = bc.turbulence_bc()
    for filename, content in turbulence_content.items():
        bc.write_bc(content, bc_path / filename)

    bc.export_bc_to_csv(bc_path.parent / "boundary_conditions_summary.csv")
    return bc
