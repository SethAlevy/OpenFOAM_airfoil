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

### Run example case
```powershell
cd /app/src/simulation/run/
bash ./run_case.sh --working-dir /app/cases --setup-file /app/examples/naca0012_study/aoa_5.json --case-name aoa_5
```

# Documentation



# Requirements

- Python 3.13
- Poetry for dependency management
- Docker

# Contact

GitHub: [Arek Drabik](https://github.com/SethAlevy)
