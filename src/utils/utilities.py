
"""
General utility functions for OpenFOAM airfoil simulations.
"""

from pathlib import Path
from io import StringIO
from typing import Optional, List
import numpy as np
import pandas as pd
import pyvista as pv
from scipy.interpolate import interp1d
from utils.logger import SimpleLogger as logger
from settings.initial_settings import InitialSettingsReader


def resolve_value(value: any, setup: dict, key: str, default: any) -> any:
    """
    Check if value is not None, return it; otherwise get from setup dict or return 
    default.

    Args:
        value (any): Value to check.
        setup (dict): Setup dictionary to get value from if value is None.
        key (str): Key to look up in setup dictionary.
        default (any): Default value to return if key not in setup.

    Returns:
        any: The resolved value.
    """
    if value is not None:
        return value
    return setup.get(key, default)


def split_naca_designation(designation: str) -> tuple[bool, str]:
    """
    Split the airfoil designation NACA type and get digits.

    Args:
        designation (str): Airfoil designation string.

    Returns:
        tuple[bool, str]: (is_naca, digits). is_naca is True if 'NACA' or 'naca' in
        designation, digits is the string of digits after 'NACA' (empty string if not
        NACA).
    """
    if "NACA" in designation or "naca" in designation:
        digits = ''.join(
            filter(str.isdigit, designation.replace("NACA", "").replace("naca", "")))
        return True, digits
    return False, ""


def read_foam_log_file(file_path: Path) -> np.ndarray:
    """
    Read a foamLog output file and return the data as numpy array.

    Args:
        file_path: Path to the log file

    Returns:
        numpy array with the data

    Raises:
        FileNotFoundError: If the log file does not exist
        ValueError: If the file cannot be parsed or is empty
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Log file not found: {file_path}")

    return np.loadtxt(file_path, comments='#')


def get_velocity_magnitude(mesh: pv.UnstructuredGrid, vel_field: str = 'U') -> np.ndarray:
    """
    Calculate velocity magnitude from velocity vector field.

    Args:
        mesh: PyVista mesh containing velocity field
        vel_field: Name of velocity field (default: 'U')

    Returns:
        Velocity magnitude array

    Raises:
        KeyError: If velocity field not found in mesh
    """
    if vel_field not in mesh.array_names:
        raise KeyError(
            f"Velocity field '{vel_field}' not found. Available: {mesh.array_names}")

    U = mesh[vel_field]
    return np.linalg.norm(U, axis=1) if U.ndim > 1 else U


def find_latest_force_coeffs_file(case_dir: Path) -> Path:
    """
    Find the most recent OpenFOAM forceCoeffs output file for a given case.

    Args:
        case_dir (Path): Path to the case directory.

    Returns:
        Path: Path to the selected coefficient.dat file.

    Raises:
        FileNotFoundError: If the forceCoeffs directory or file cannot be found.
    """
    force_coeffs_dir = case_dir / "postProcessing" / "forceCoeffs"
    candidates = list(force_coeffs_dir.rglob("coefficient.dat"))
    if not candidates:
        raise FileNotFoundError(f"No coefficient.dat found under: {force_coeffs_dir}")

    times = []
    for p in candidates:
        t = float(p.parent.name)
        times.append((t, p.stat().st_mtime, p))
    times.sort(reverse=True)
    return times[0][2]


def read_force_coeffs_dat(file_path: Path) -> dict[str, np.ndarray]:
    """
    Read an OpenFOAM coefficient.dat file into named columns.

    Args:
        file_path (Path): Path to coefficient.dat.

    Returns:
        dict[str, np.ndarray]: Mapping from column name to data array.

    Raises:
        ValueError: If the file cannot be parsed into numeric data.
    """
    header_cols: Optional[List[str]] = None
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s.startswith("#"):
                continue

            maybe = s.lstrip("#").strip()
            parts = maybe.split()
            if parts and parts[0] == "Time":
                header_cols = parts

    data = np.loadtxt(file_path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if header_cols is None:
        header_cols = ["Time", "Cd", "Cs", "Cl", "CmRoll", "CmPitch", "CmYaw"]

    ncols = min(len(header_cols), data.shape[1])
    return {header_cols[i]: data[:, i] for i in range(ncols)}


def combine_postprocessing_summaries(
        working_path: Path,
        output_filename: str = "combined_postprocessing_summary.csv"
) -> None:
    """
    Combine all postprocessing_summary.csv files from subdirectories of working_path
    into a single CSV, sorted by AngleOfAttack.

    Args:
        working_path (Path): Path containing case directories with
            postprocessing_summary.csv files.
        output_filename (str): Name of the combined CSV file to create.
    """
    working_path = Path(working_path)
    summary_files = list(working_path.rglob("postprocessing_summary.csv"))
    if not summary_files:
        logger.log(f"No postprocessing_summary.csv files found in {working_path}")
        return

    dfs = []
    for file in summary_files:
        try:
            df = pd.read_csv(file)
            df["CaseDir"] = str(file.parent)
            dfs.append(df)
        except Exception as e:
            logger.log(f"Failed to read {file}: {e}")

    if not dfs:
        logger.log("No valid summary files to combine.")
        return

    combined_df = pd.concat(dfs, ignore_index=True)
    if "AngleOfAttack" in combined_df.columns:
        combined_df = combined_df.sort_values("AngleOfAttack")
    output_path = working_path / output_filename
    combined_df.to_csv(output_path, index=False)
    logger.log(f"Combined summary saved to {output_path}")


def export_postprocessing_summary(
    case_dir: Path,
    output_csv: Path = None,
    mean_n: int = 50
) -> None:
    """
    Export postprocessing summary CSV for an OpenFOAM airfoil case.

    Args:
        case_dir (Path): Path to the case directory containing the results.
        output_csv (Path): Output CSV file path (default:
                <case_dir>/postprocessing_summary.csv).
        mean_n (int): Number of last timesteps to average for coefficients (default: 50).
    """
    output_csv = output_csv or (case_dir / "postprocessing_summary.csv")

    json_file = case_dir / f"{case_dir.name}.json"
    settings_reader = InitialSettingsReader(json_file)

    bc_csv = case_dir / "boundary_conditions_summary.csv"
    coeff_file = find_latest_force_coeffs_file(case_dir)
    coeffs = read_force_coeffs_dat(coeff_file)

    designation = settings_reader.airfoil_settings.get("Designation", "")
    aoa = settings_reader.airfoil_settings.get("AngleOfAttack", "")

    df_bc = pd.read_csv(bc_csv)
    bc = dict(zip(df_bc['Parameter'], df_bc['Value']))
    mach = bc.get("Mach Number", "")
    reynolds = bc.get("Reynolds Number", "")
    velocity = bc.get("Velocity", "")
    density = bc.get("Density", "")
    chord = bc.get("Chord Length", "")

    n = len(coeffs["Time"])
    idx_start = max(0, n - mean_n)
    cd = np.mean(coeffs["Cd"][idx_start:])
    cl = np.mean(coeffs["Cl"][idx_start:])
    cmpitch = np.mean(coeffs.get("CmPitch", coeffs.get("Cm", np.zeros(n)))[idx_start:])

    df = pd.DataFrame([{
        "Designation": designation,
        "AngleOfAttack": aoa,
        "Cd": cd,
        "Cl": cl,
        "CmPitch": cmpitch,
        "Mach": mach,
        "Reynolds": reynolds,
        "Velocity": velocity,
        "Density": density,
        "Chord": chord
    }])

    df.to_csv(output_csv, index=False)
    logger.log(f"Exported summary to {output_csv}")


def read_uiuc_format_reference_csv(csv_path : Path) -> pd.DataFrame:
    """
    Reads a UIUC/Xfoil-style reference CSV, skipping header lines. This is the 
    default format used inside the repository.

    Args:
        csv_path (Path): Path to the reference CSV file.

    Returns:
        pd.DataFrame: DataFrame containing the airfoil data.

    """
    with open(csv_path, 'r') as f:
        lines = f.readlines()
    # Find the line where the data table starts
    for idx, line in enumerate(lines):
        if line.strip().startswith("Alpha,"):
            start_idx = idx
            break
    else:
        raise ValueError(f"Could not find data table in {csv_path}")
    data_str = "".join(lines[start_idx:])
    return pd.read_csv(StringIO(data_str))


def interpolate_airfoil_coefficients(
        csv_path: Path,
        angle_of_attack: float,
        min_angle: float = -5.0,
        max_angle: float = 17.0
) -> dict[str, float]:
    """
    Interpolate aerodynamic coefficients (Cl, Cd, Cm) from a UIUC/Xfoil-style CSV
    for a given angle of attack in the range.

    Args:
        csv_path (Path): Path to the reference CSV file.
        angle_of_attack (float): Angle of attack in degrees.
        min_angle (float): Minimum valid angle of attack (default: -5.0).
        max_angle (float): Maximum valid angle of attack (default: 17.0).

    Returns:
        dict[str, float]: Dictionary containing interpolated 'Cl', 'Cd', and 'Cm' values.

    Raises:
        ValueError: If angle of attack is outside the valid range (-5 to +17 degrees).
        FileNotFoundError: If CSV file not found.
        KeyError: If required columns are missing in CSV.
    """
    if not isinstance(csv_path, Path):
        csv_path = Path(csv_path)

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    if angle_of_attack < min_angle or angle_of_attack > max_angle:
        raise ValueError(
            f"Angle of attack {angle_of_attack}° is outside valid range "
            f"[{min_angle}°, {max_angle}°]"
        )

    df = read_uiuc_format_reference_csv(csv_path)

    required_cols = ['Alpha', 'Cl', 'Cd', 'Cm']
    if missing_cols := [col for col in required_cols if col not in df.columns]:
        raise KeyError(f"Missing columns in CSV: {missing_cols}")

    alpha = df['Alpha'].values
    cl = df['Cl'].values
    cd = df['Cd'].values
    cm = df['Cm'].values

    interp_cl = interp1d(alpha, cl, kind='cubic', fill_value='extrapolate')
    interp_cd = interp1d(alpha, cd, kind='cubic', fill_value='extrapolate')
    interp_cm = interp1d(alpha, cm, kind='cubic', fill_value='extrapolate')

    cl_interp = float(interp_cl(angle_of_attack))
    cd_interp = float(interp_cd(angle_of_attack))
    cm_interp = float(interp_cm(angle_of_attack))

    return {
        'Cl': cl_interp,
        'Cd': cd_interp,
        'Cm': cm_interp,
        'Alpha': angle_of_attack
    }
