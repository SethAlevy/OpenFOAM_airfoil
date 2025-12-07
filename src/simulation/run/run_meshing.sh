#!/bin/bash

setup_file="$1"

set --

source /usr/lib/openfoam/openfoam2406/etc/bashrc || true
set -e

if [ -z "$setup_file" ]; then
    echo "Usage: $0 <setup_file>"
    exit 1
fi

mesher=$(cat "$setup_file" | tr -d '\r' | jq -r '.Mesh.Mesher // "cfMesh"')
n_procs=$(cat "$setup_file" | tr -d '\r' | jq -r '.Simulation.Decomposition.NumberOfSubdomains // 1')
print_logs=$(cat "$setup_file" | tr -d '\r' | jq -r '.Simulation.PrintLogs // false')

touch open.foam
echo "Selected mesher: ${mesher}"
echo "Working dir: $PWD"
echo "Print logs: ${print_logs}"

mkdir -p mesh_logs

run_cmd() {
    local cmd="$1"
    local logfile="mesh_logs/$2"
    
    if [ "$print_logs" = "true" ]; then
        $cmd 2>&1 | tee "$logfile"
    else
        $cmd > "$logfile" 2>&1
    fi
}

if [ "${mesher}" = "cfMesh" ]; then
    
    # This workflow uses a single call to cartesianMesh.
    # The python script prepares a composite STL file defining the whole domain.
    echo "Running cfMesh (cartesian2DMesh) with composite STL..."
    run_cmd "cartesian2DMesh" "cartesian2DMesh.log"

    # The extrudeMesh step is not needed as the composite STL is already 3D.

else
    echo "Running SnappyHexMesh pipeline ..."
    if [ -f "system/surfaceFeaturesDict" ]; then
        echo "surfaceFeaturesDict found. Running surfaceFeatures ..."
        run_cmd "surfaceFeatures" "surfaceFeatures.log"
    fi

    echo "Running blockMesh ..."
    run_cmd "blockMesh" "blockMesh.log"

    if [ "$n_procs" -le 1 ]; then
        echo "Single processor. Running snappyHexMesh ..."
        run_cmd "snappyHexMesh -overwrite" "snappyHexMesh.log"
    else
        echo "Decomposing for parallel snappyHexMesh ..."
        run_cmd "decomposePar -force" "decomposePar.log"
        echo "Running snappyHexMesh in parallel ..."
        run_cmd "mpirun -np $n_procs snappyHexMesh -parallel -overwrite" "snappyHexMesh.log"
        echo "Reconstructing ..."
        run_cmd "reconstructPar -constant" "reconstructPar.log"
        rm -rf processor*
    fi
fi

echo "Renumbering ..."
run_cmd "renumberMesh -overwrite" "renumberMesh.log"

echo "Checking mesh ..."
run_cmd "checkMesh" "checkMesh.log"

echo "Meshing completed. Logs saved in mesh_logs/"
