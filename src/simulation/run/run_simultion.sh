#!/bin/bash

source /opt/openfoam11/etc/bashrc

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <n_processors>"
    exit 1
fi

n_procs="$1"

solver=$(grep "^application" system/controlDict | awk '{print $2}')

if [ -z "$solver" ]; then
    echo "Solver not found in system/controlDict. Using default: simpleFoam"
    exit 1
fi

echo "Starting the simulation with $n_procs processors in $PWD"
echo "Decomposing the domain for parallel simulation..."
decomposePar -force 2>&1 | tee decomposePar.log

echo "Running the solver ($solver) in parallel..."
mpirun -np "$n_procs" "$solver" -parallel 2>&1 | tee "$solver".log

echo "Reconstructing the results..."
reconstructPar 2>&1 | tee reconstructPar.log
rm -rf processor*

foamLog "$solver".log

mkdir simulation_logs
mv *.log simulation_logs/

echo "Simulation completed. Logs saved in simulation_logs/."
