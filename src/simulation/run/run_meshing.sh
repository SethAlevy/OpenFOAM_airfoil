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

touch open.foam
echo "Selected mesher: ${mesher}"
echo "Working dir: $PWD"

if [ "${mesher}" = "cfMesh" ]; then
    
    # This workflow uses a single call to cartesianMesh.
    # The python script prepares a composite STL file defining the whole domain.
    echo "Running cfMesh (cartesianMesh) with composite STL..."
    cartesian2DMesh 2>&1 | tee cartesian2DMesh.log

    # The extrudeMesh step is not needed as the composite STL is already 3D.

else
    echo "Running SnappyHexMesh pipeline ..."
    if [ -f "system/surfaceFeaturesDict" ]; then
        echo "surfaceFeaturesDict found. Running surfaceFeatures ..."
        surfaceFeatures 2>&1 | tee surfaceFeatures.log
    fi

    echo "Running blockMesh ..."
    blockMesh 2>&1 | tee blockMesh.log

    if [ "$n_procs" -le 1 ]; then
        echo "Single processor. Running snappyHexMesh ..."
        snappyHexMesh -overwrite 2>&1 | tee snappyHexMesh.log
    else
        echo "Decomposing for parallel snappyHexMesh ..."
        decomposePar -force 2>&1 | tee decomposePar.log
        echo "Running snappyHexMesh in parallel ..."
        mpirun -np "$n_procs" snappyHexMesh -parallel -overwrite 2>&1 | tee snappyHexMesh.log
        echo "Reconstructing ..."
        reconstructPar -constant 2>&1 | tee reconstructPar.log
        rm -rf processor*
    fi
fi

echo "Renumbering ..."
renumberMesh -overwrite 2>&1 | tee renumberMesh.log

echo "Checking mesh ..."
checkMesh 2>&1 | tee checkMesh.log

mkdir -p mesh_logs
mv *.log mesh_logs/
echo "SnappyHexMesh meshing completed."
