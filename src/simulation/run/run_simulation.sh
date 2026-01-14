#!/bin/bash

setup_file="$1"

set --

source /usr/lib/openfoam/openfoam2406/etc/bashrc || true
set -e

if [ -z "$setup_file" ]; then
    echo "Usage: $0 <setup_file>"
    exit 1
fi

# Extract solver from system/controlDict
if [ -f "system/controlDict" ]; then
    solver=$(foamDictionary system/controlDict -entry application -value 2>/dev/null || echo "simpleFoam")
else
    solver="simpleFoam"
fi

n_procs=$(cat "$setup_file" | tr -d '\r' | jq -r '.Simulation.Decomposition.NumberOfSubdomains // 1')
print_logs=$(cat "$setup_file" | tr -d '\r' | jq -r '.Simulation.PrintLogs // false')

echo "Number of processors: ${n_procs}"
echo "Print logs: ${print_logs}"

mkdir -p simulation_logs

run_cmd() {
    local cmd="$1"
    local logfile="simulation_logs/$2"
    
    if [ "$print_logs" == "true" ]; then
        $cmd 2>&1 | tee "$logfile"
    else
        $cmd > "$logfile" 2>&1
    fi
}

if [ "$n_procs" -le 1 ]; then
    echo "Single processor detected. Running $solver without decomposition..."
    run_cmd "$solver" "${solver}.log"
else
    echo "Starting the simulation with $n_procs processors in $PWD"
    echo "Decomposing the domain for parallel simulation..."
    run_cmd "decomposePar -force" "decomposePar.log"

    echo "Running the solver $solver in parallel..."
    run_cmd "mpirun -np $n_procs $solver -parallel" "${solver}.log"

    echo "Reconstructing the results..."
    run_cmd "reconstructPar" "reconstructPar.log"
    rm -rf processor*
fi

echo "Extracting logs..."
foamLog simulation_logs/${solver}.log

echo "Simulation completed. Logs saved in simulation_logs/."
