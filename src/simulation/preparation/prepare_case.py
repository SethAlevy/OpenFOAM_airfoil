from pathlib import Path
from simulation.openfoam.mesh import block_mesh_dict, snappy_hex_mesh_dict
from simulation.openfoam.system_dir import control_dict, fv_solution_dict, \
      fv_schemes_dict
from simulation.openfoam.constant_dir import transport_properties_dict, \
      turbulence_properties_dict
from simulation.openfoam.boundary_condition import BoundaryConditions
from templates.airfoil_template import Airfoil
from templates.initial_settings_template import Settings


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
    create_meshing_files(working_path / case_name, airfoil, setup, mesher)
    create_system_files(working_path / case_name, setup)
    create_constant_files(working_path / case_name, setup)
    create_boundary_conditions_files(working_path / case_name / "0", setup)

    airfoil.to_stl(working_path / case_name / "constant" / "triSurface" / "airfoil.stl")


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
        mesher (str): The meshing method to use (e.g., "blockMesh", "snappyHexMesh").
    """

    block_mesh_dict(airfoil, setup, (case_path / "system"))

    if mesher == "snappyHexMesh":
        snappy_hex_mesh_dict(airfoil, setup, (case_path / "system"))


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


def create_constant_files(case_path: Path, setup: Settings) -> None:
    """
    Creates constant directory files for the OpenFOAM case.

    Args:
        case_path (Path): The path to the case directory.
        setup (Settings): The simulation settings.
    """
    transport_properties_dict(setup, (case_path / "constant" / "transportProperties"))
    turbulence_properties_dict(setup, (case_path / "constant" / "turbulenceProperties"))


def create_boundary_conditions_files(case_path: Path, setup: Settings) -> None:
    """
    Creates boundary condition files for the OpenFOAM case.

    Args:
        case_path (Path): The path to the case directory.
        setup (Settings): The simulation settings.
    """
    bc = BoundaryConditions(setup=setup)
    velocity_content = bc.velocity_bc()
    pressure_content = bc.pressure_bc()

    bc.write_bc(velocity_content, case_path / "U")
    bc.write_bc(pressure_content, case_path / "p")
