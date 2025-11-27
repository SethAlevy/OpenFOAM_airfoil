#!/bin/bash

source /opt/openfoam11/etc/bashrc

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <n_processors>"
    exit 1
fi

n_procs="$1"

echo "Starting the meshing process with $n_procs processors in $PWD"

if [ -f "system/surfaceFeatureExtractDict" ]; then
    echo "surfaceFeatureExtractDict found. Running surfaceFeatureExtract..."
    surfaceFeatureExtract 2>&1 | tee surfaceFeatureExtract.log
fi

echo "Running blockMesh..."
blockMesh 2>&1 | tee blockMesh.log

echo "Decomposing the domain for parallel snappyHexMesh..."
decomposePar -force 2>&1 | tee decomposePar.log

echo "Running snappyHexMesh in parallel..."
mpirun -np "$n_procs" snappyHexMesh -parallel -overwrite 2>&1 | tee snappyHexMesh.log

echo "Reconstructing the mesh..."
reconstructParMesh -constant 2>&1 | tee reconstructParMesh.log
rm -rf processor*

echo "Renumbering the mesh..."
renumberMesh -overwrite 2>&1 | tee renumberMesh.log

echo "Checking the mesh quality..."
checkMesh 2>&1 | tee checkMesh.log

mkdir mesh_logs
mv *.log mesh_logs/

echo "Meshing process completed. Logs saved in mesh_logs/."