#!/bin/bash

setup_file="$1"

set --

source /usr/lib/openfoam/openfoam2406/etc/bashrc || true
set -ex

if [ -z "$setup_file" ]; then
    echo "Usage: $0 <setup_file>"
    exit 1
fi

# Check if boundary_conditions.csv exists
bc_csv="boundary_conditions_summary.csv"
if [ ! -f "$bc_csv" ]; then
    echo "Error: boundary_conditions_summary.csv not found in case directory"
    exit 1
fi

# Extract solver from system/controlDict
if [ -f "system/controlDict" ]; then
    solver=$(foamDictionary system/controlDict -entry application -value 2>/dev/null || echo "simpleFoam")
else
    solver="simpleFoam"
fi

# Extract parameters from CSV file
chord=$(cat "$setup_file" | tr -d '\r' | jq -r '.Airfoil.Chord // 1.0')
density=$(grep "^Density," "$bc_csv" | cut -d',' -f3)
velocity=$(grep "^Velocity," "$bc_csv" | cut -d',' -f3)

# Extract angle of attack from JSON (not in CSV)
aoa=$(cat "$setup_file" | tr -d '\r' | jq -r '.Airfoil.AngleOfAttack // 0.0')

# For 2D simulation, use unit span
z_span=1.0
airfoil_patch="airfoil"

# Calculate reference area (chord * unit span for 2D)
Aref=$(echo "$chord * $z_span" | bc -l)

# Convert angle of attack to radians for lift/drag direction
aoa_rad=$(echo "$aoa * 3.14159265359 / 180" | bc -l)
lift_x=$(echo "-1 * s($aoa_rad)" | bc -l)
lift_y=$(echo "c($aoa_rad)" | bc -l)
drag_x=$(echo "c($aoa_rad)" | bc -l)
drag_y=$(echo "s($aoa_rad)" | bc -l)
cofr_x=$(echo "$chord * 0.25" | bc -l)

echo "Post-processing with parameters:"
echo "  Solver: $solver"
echo "  Chord: $chord m"
echo "  Span (2D): $z_span m"
echo "  Reference Area: $Aref m²"
echo "  Velocity: $velocity m/s"
echo "  Density: $density kg/m³"
echo "  Angle of Attack: $aoa degrees"
echo "  Airfoil Patch: $airfoil_patch"

echo "Calculating force coefficients..."
postProcess -time "0:" -func "forceCoeffs(libs=(forces), patches=($airfoil_patch), rho=rhoInf, rhoInf=$density, CofR=($cofr_x 0 0), liftDir=($lift_x $lift_y 0), dragDir=($drag_x $drag_y 0), pitchAxis=(0 0 1), magUInf=$velocity, lRef=$chord, Aref=$Aref)"

echo "Calculating forces..."
postProcess -latestTime -func "forces(libs=(forces), patches=($airfoil_patch), rho=rhoInf, rhoInf=$density, CofR=($cofr_x 0 0))"

echo "Calculating pressure coefficient..."
postProcess -latestTime -func "pressureCoefficient(U=U, rho=rhoInf, rhoInf=$density, pInf=0, UInf=($velocity 0 0))"

echo "Calculating y+ values..."
$solver -postProcess -latestTime -func yPlus

echo "Calculating wall shear stress..."
$solver -postProcess -latestTime -func wallShearStress

echo "Converting to VTK..."
foamToVTK -latestTime

echo "Generating visualization plots..."
poetry run python3 /app/src/postprocess/generate_plots.py \
    --case-dir "$(pwd)" \
    --plots all \
    --formats png html

echo "All post-processing and visualization completed."
