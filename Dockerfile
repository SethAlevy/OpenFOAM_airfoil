# Dockerfile for OpenFOAM Airfoil Project
# Ubuntu base image with OpenFOAM installation

FROM openfoam/openfoam11-paraview510:latest

# The official OpenFOAM 11 image with ParaView 5.10 included.
# For more info: https://hub.docker.com/r/openfoam/openfoam11-paraview510

# Default command
CMD ["/bin/bash"]
