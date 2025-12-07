#!/bin/bash

setup_file="$1"

set --

source /usr/lib/openfoam/openfoam2406/etc/bashrc || true
set -e

if [ -z "$setup_file" ]; then
    echo "Usage: $0 <setup_file>"
    exit 1
fi

n_procs=$(cat "$setup_file" | tr -d '\r' | jq -r '.Simulation.Decomposition.NumberOfSubdomains // 1')
print_logs=$(cat "$setup_file" | tr -d '\r' | jq -r '.Simulation.PrintLogs // true')

echo "Number of processors: ${n_procs}"
echo "Print logs: ${print_logs}"

# Create log directory
mkdir -p simulation_logs

# Function to run command with or without printing
run_cmd() {
    local cmd="$1"
    local logfile="simulation_logs/$2"
    
    if [ "$print_logs" = "true" ]; then
        $cmd 2>&1 | tee "$logfile"
    else
        $cmd > "$logfile" 2>&1
    fi
}

if [ "$n_procs" -le 1 ]; then
    echo "Single processor detected. Running simulation without decomposition..."
    run_cmd "simpleFoam" "simpleFoam.log"
else
    echo "Starting the simulation with $n_procs processors in $PWD"
    echo "Decomposing the domain for parallel simulation..."
    run_cmd "decomposePar -force" "decomposePar.log"

    echo "Running the solver (simpleFoam) in parallel..."
    run_cmd "mpirun -np $n_procs simpleFoam -parallel" "simpleFoam.log"

    echo "Reconstructing the results..."
    run_cmd "reconstructPar" "reconstructPar.log"
    rm -rf processor*
fi

echo "Extracting logs..."
foamLog simulation_logs/simpleFoam.log

echo "Simulation completed. Logs saved in simulation_logs/."
