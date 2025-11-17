# OpenFOAM Airfoil Simulation Project

This project provides a Docker environment for generating and simulating airfoils using OpenFOAM on Ubuntu.

## Getting Started

### Build the Docker Image
```powershell
docker build -t openfoam-airfoil .
```

### Run the Container
```powershell
docker run -it openfoam-airfoil
```

This will drop you into a bash shell with OpenFOAM sourced and ready to use.

## Files
- `Dockerfile`: Builds the Ubuntu + OpenFOAM environment.
- `.dockerignore`: Prevents unnecessary files from being added to the image.
- `README.md`: Project instructions.

## Next Steps
- Add scripts for airfoil geometry generation.
- Automate OpenFOAM case setup and simulation.

---

For more information on OpenFOAM, visit: https://openfoam.org/
