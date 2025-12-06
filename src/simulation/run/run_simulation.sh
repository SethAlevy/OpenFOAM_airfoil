#!/bin/bash

setup_file="$1"

set --

source /usr/lib/openfoam/openfoam2406/etc/bashrc || true
set -ex

if [ -z "$setup_file" ]; then
    echo "Usage: $0 <setup_file>"
    exit 1
fi

n_procs=$(cat "$setup_file" | tr -d '\r' | jq -r '.Simulation.Decomposition.NumberOfSubdomains // 1')

if [ "$n_procs" -le 1 ]; then
    echo "Single processor detected. Running simulation without decomposition..."
    simpleFoam 2>&1 | tee simpleFoam.log
else
    echo "Starting the simulation with $n_procs processors in $PWD"
    echo "Decomposing the domain for parallel simulation..."
    decomposePar -force 2>&1 | tee decomposePar.log

    echo "Running the solver (simpleFoam) in parallel..."
    mpirun -np "$n_procs" simpleFoam -parallel 2>&1 | tee simpleFoam.log

    echo "Reconstructing the results..."
    reconstructPar 2>&1 | tee reconstructPar.log
    # rm -rf processor*
fi

foamLog simpleFoam.log

mkdir -p simulation_logs
mv *.log simulation_logs/

echo "Simulation completed. Logs saved in simulation_logs/."
