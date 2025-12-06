import argparse
from pathlib import Path
from settings.initial_settings import InitialSettingsReader
import airfoil.airfoil as af
from simulation.preparation.case_structure import prepare_openfoam_case
import utils.utilities as ut
from utils.logger import SimpleLogger


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Prepare OpenFOAM case structure and files based on the passed"
        " configuration."
    )

    parser.add_argument(
        '--working-dir',
        type=Path,
        required=False,
        default=Path(__file__).parent.parent.parent.parent.resolve(),
        help='Path to the directory where the case will be prepared. Main repository'
        ' path is the default.'
    )

    parser.add_argument(
        '--setup-file',
        type=Path,
        required=True,
        help='Path to the setup configuration json file. The file contains all'
        ' parameters that define the simulation.'
    )

    parser.add_argument(
        '--case-name',
        type=str,
        required=False,
        default=None,
        help='Name of the case directory to be created within the working directory.'
    )

    return parser.parse_args()


def prepare_case():
    args = parse_arguments()
    working_dir = args.working_dir
    setup_file = args.setup_file

    setup = InitialSettingsReader(setup_file)

    airfoil_settings = setup.airfoil_settings

    designation = airfoil_settings.get("Designation")

    case_name = args.case_name if args.case_name is not None else f"case_{designation}"

    if bool(airfoil_settings.get("GenerateNACA")):
        naca, digits = ut.split_naca_designation(designation)
        if not naca:
            raise ValueError(
                f"Invalid NACA designation '{designation}' for airfoil generation."
            )
        
        if len(digits) == 4:
            airfoil = af.NACA4(designation=digits, setup=setup)
        elif len(digits) == 5:
            airfoil = af.NACA5(designation=digits, setup=setup)
        else:
            raise ValueError(
                "Unsupported number of digits. Only NACA 4 and 5 digit airfoils"
                " available."
            )

        SimpleLogger.log(f"Generated airfoil: {airfoil}")
    
    elif bool(airfoil_settings.get("DownloadUIUC")):
        airfoil = af.UIUCAirfoil(designation=designation, setup=setup)

        SimpleLogger.log(f"Downloaded airfoil: {airfoil}")

    elif bool(airfoil_settings.get("LoadFromFile")):
        file_path = airfoil_settings.get("LoadFilePath", None)
        if file_path is None or not Path(file_path).is_file():
            raise FileNotFoundError(
                f"Airfoil file path '{file_path}' is invalid or does not exist."
            )
        # TODO dat file loader

        SimpleLogger.log(f"Loaded airfoil from file: {airfoil}")

    SimpleLogger.log(f"Preparing case '{case_name}' in '{working_dir}'")
    prepare_openfoam_case(
        working_path=working_dir,
        case_name=case_name,
        airfoil=airfoil,
        setup=setup
    )


if __name__ == "__main__":
    prepare_case()
