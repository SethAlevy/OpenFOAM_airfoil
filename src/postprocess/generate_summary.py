import argparse
import json
import pandas as pd
from pathlib import Path
import numpy as np
from utils.logger import SimpleLogger as logger
from utils.utilities import find_latest_force_coeffs_file, read_force_coeffs_dat


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Export postprocessing summary CSV for an OpenFOAM airfoil case."
    )
    parser.add_argument(
        '--case-dir',
        type=Path,
        required=True,
        help='Path to the case directory containing the results.'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output CSV file path (default: <case_dir>/postprocessing_summary.csv)'
    )
    parser.add_argument(
        '--mean-n',
        type=int,
        default=1,
        help='Number of last timesteps to average for coefficients (default: 1, i.e. last value only).'
    )
    return parser.parse_args()


def read_bc_csv(bc_csv: Path):
    df = pd.read_csv(bc_csv)
    bc = dict(zip(df['Parameter'], df['Value']))
    mach = bc.get("Mach Number", "")
    reynolds = bc.get("Reynolds Number", "")
    velocity = bc.get("Velocity", "")
    density = bc.get("Density", "")
    chord = bc.get("Chord Length", "")
    return mach, reynolds, velocity, density, chord


def generate_summary():
    args = parse_arguments()
    case_dir = args.case_dir
    output_csv = args.output or (case_dir / "postprocessing_summary.csv")
    mean_n = args.mean_n

    json_file = case_dir / f"{case_dir.name}.json"
    bc_csv = case_dir / "boundary_conditions_summary.csv"
    coeff_file = find_latest_force_coeffs_file(case_dir)
    coeffs = read_force_coeffs_dat(coeff_file)

    with open(json_file, "r") as f:
        j = json.load(f)
    designation = j["Airfoil"].get("Designation", "")
    aoa = j["Airfoil"].get("AngleOfAttack", "")

    mach, reynolds, velocity, density, chord = read_bc_csv(bc_csv)

    # Get mean or last value for coefficients
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


if __name__ == "__main__":
    generate_summary()
