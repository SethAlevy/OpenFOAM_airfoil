import json
from utils.logger import SimpleLogger
from pathlib import Path


class InitialSettingsReader():
    def __init__(self, filepath: Path):
        """
        Class to read and validate the initial settings json file.

        Args:
            filepath (Path): Path to the file.
        """
        self.filepath = filepath
        self.settings = self.read_settings()
        self.validate_settings()

    def read_settings(self) -> dict:
        """
        Read the settings json file.

        Returns:
            dict: settings json file parsed as dictionary.
        """
        with open(self.filepath, 'r') as file:
            return json.load(file)

    @property
    def airfoil_settings(self) -> dict:
        """
        Contains parameters related to airfoil generation, like designation, chord.

        Returns:
            dict: Airfoil parameters.
        """
        return self.settings.get("Airfoil", {})

    @property
    def mesh_settings(self) -> dict:
        """
        Contains parameters related to mesh generation in openFOAM. It allows to handle
        the creation of blockMesh and further refinement using snappyHexMesh or cfMesh.

        Returns:
            dict: Mesh parameters.
        """
        return self.settings.get("Mesh", {})

    @property
    def simulation_settings(self) -> dict:
        """
        Contains parameters to handle the physical properties of the simulation 
        like turbulence model, fluid properties and openFOAM files like controlDict,
        fvSolution, fvSchemes and boundary conditions.

        Returns:
            dict: Simulation parameters.
        """
        return self.settings.get("Simulation", {})

    @staticmethod
    def recursive_check(template: dict, actual: dict, path: str = ""):
        """
        Performs recursive check on the actual json and template file. Check if required
        keys are present, the keys are in the right position and warns for extra keys
        in case of mismatch.

        Args:
            template (dict): Template dictionary to check against.
            actual (dict): Actual dictionary to validate.
            path (str, optional): Current path in the dictionary for nested keys.
        """
        # Warn for missing high-level keys
        if path == "":
            for key in template:
                if key not in actual:
                    raise KeyError(f"Missing top-level key '{key}' in settings.")
        # Check for misplaced or extra keys
        for key in actual:
            if key not in template:
                SimpleLogger.warning(
                    f"Key '{path + key}' is not allowed in this position.")
            elif isinstance(actual[key], dict) and isinstance(template[key], dict):
                InitialSettingsReader.recursive_check(
                    template[key], actual[key], path + key + ".")
            elif isinstance(actual[key], dict) and not isinstance(template[key], dict):
                SimpleLogger.warning(
                    f"Key '{path + key}' should not be a dictionary here.")

    def validate_settings(self) -> None:
        """
        Validate the settings against the template structure. Logs warnings for wrong
        placed or extra keys.
        """
        template = self.template_settings()
        self.recursive_check(template, self.settings)

    @staticmethod
    def template_settings() -> dict:
        """
        A full settings template to allow validation checks on the settings json.

        Returns:
            dict: Template settings dictionary.
        """
        return {
            "Airfoil": {
                "Designation": "",
                "GenerateNACA4": False,
                "GenerateNACA5": False,
                "DownloadUIUC": False,
                "LoadFromFile": False,
                "LoadFilePath": "",
                "Chord": 1.0,
                "Resolution": 200
            },
            "Mesh": {
                "Mesher": "",
                "BlockMesh": {
                    "BoundingBox": {
                        "InletDistance": 0.0,
                        "OutletDistance": 0.0,
                        "TopDistance": 0.0,
                        "BottomDistance": 0.0,
                        "BaseCellSize": 0.0,
                        "ZMin": 0.0,
                        "ZMax": 0.0,
                        "Nz": 1
                    }
                },
                "SnappyHexMesh": {
                    "Geometry": {
                        "RefinementBoxLevel": 0,
                        "SphereTipLevel": 0,
                        "SphereTipRadius": 0.0,
                        "SphereLeadingEdgeLevel": 0,
                        "SphereLeadingEdgeRadius": 0.0
                    },
                    "SnapControls": {
                        "NSmoothPatch": 0,
                        "Tolerance": 0.0,
                        "NSolveIter": 0,
                        "NRelaxIter": 0,
                        "NFeatureSnapIter": 0
                    },
                    "CastellatedMeshControls": {
                        "MaxLocalCells": 0,
                        "MaxGlobalCells": 0,
                        "MinRefinementCells": 0,
                        "MaxLoadUnbalance": 0.0,
                        "NCellsBetweenLevels": 0,
                        "MinSurfaceRefinementLevel": 0,
                        "MaxSurfaceRefinementLevel": 0,
                        "FeatureRefinementLevel": 0,
                        "ResolveFeatureAngle": 0,
                        "MakeBaffle": False,
                        "AllowFreeStandingZoneFaces": False
                    },
                    "AddLayersControls": {
                        "AddLayers": False,
                        "NSurfaceLayers": 0,
                        "ExpansionRatio": 0.0,
                        "FinalLayerThickness": 0.0,
                        "MinThickness": 0.0,
                        "NGrow": 0,
                        "FeatureAngle": 0,
                        "NRelaxIter": 0,
                        "NSmoothSurfaceNormals": 0,
                        "NSmoothNormals": 0,
                        "NSmoothThickness": 0,
                        "MaxFaceThicknessRatio": 0.0,
                        "MaxThicknessToMedialRatio": 0.0,
                        "MinMedianAxisAngle": 0,
                        "NBufferCellsNoExtrude": 0,
                        "NLayerIter": 0
                    }
                },
                "CfMesh": {
                    "BaseCellSize": 0.0,
                    "DomainMin": 0.0,
                    "DomainMax": 0.0,
                    "AirfoilRefinementMinLevel": 0,
                    "AirfoilRefinementMaxLevel": 0,
                    "ResolveFeatureAngle": 0,
                    "CurvatureRefinementLevel": 0,
                    "ExtrudeThickness": 0.0,
                    "ExtrudeLayers": 0
                }
            },
            "Simulation": {
                "CaseName": "",
                "ControlDict": {
                    "Solver": "",
                    "StartFrom": "",
                    "StartTime": 0,
                    "EndTime": 0,
                    "DeltaT": 0,
                    "WriteControl": "",
                    "WriteInterval": 0,
                    "PurgingInterval": 0
                },
                "FvSolution": {
                    "k": {
                        "Solver": "",
                        "Smoother": "",
                        "Tolerance": 0.0,
                        "RelTol": 0.0
                    },
                    "Omega": {
                        "Solver": "",
                        "Smoother": "",
                        "Tolerance": 0.0,
                        "RelTol": 0.0
                    }
                },
                "FvSchemes": {
                    "DivSchemes": {
                        "Default": "",
                        "DivPhiU": "",
                        "DivPhiK": "",
                        "DivPhiEpsilon": "",
                        "DivPhiOmega": ""
                    },
                    "LaplacianSchemes": {
                        "Default": ""
                    },
                    "InterpolationSchemes": {
                        "Default": ""
                    }
                },
                "Decomposition": {
                    "Method": "",
                    "NumberOfSubdomains": 0
                },
                "SurfaceFeature": {
                    "IncludedAngle": 0
                },
                "TurbulenceProperties": {
                    "Turbulence": "",
                    "Model": "",
                    "Kappa": 0.0,
                    "E": 0.0,
                    "TurbulenceIntensity": 0.0,
                    "TurbulenceLengthScale": 0.0,
                    "kOmegaSST": {
                        "betaStar": 0.0,
                        "a1": 0.0,
                        "c1": 0.0,
                        "c2": 0.0,
                        "alphaK1": 0.0,
                        "alphaK2": 0.0,
                        "alphaOmega1": 0.0,
                        "alphaOmega2": 0.0,
                        "beta1": 0.0,
                        "beta2": 0.0,
                        "gamma1": 0.0,
                        "gamma2": 0.0
                    },
                    "kEpsilon": {
                        "Cmu": 0.0,
                        "C1Epsilon": 0.0,
                        "C2Epsilon": 0.0,
                        "SigmaK": 0.0,
                        "SigmaEpsilon": 0.0
                    }
                },
                "Fluid": {
                    "KinematicViscosity": 0.0,
                    "Density": 0.0,
                    "Temperature": 0.0
                },
                "BoundaryConditions": {
                    "Velocity": 0.0,
                    "Pressure": 0.0,
                    "MachNumber": 0.0,
                    "ReynoldsNumber": 0.0,
                    "InletPatchName": "",
                    "OutletPatchName": "",
                    "LowerWallPatchName": "",
                    "UpperWallPatchName": "",
                    "FrontPatchName": "",
                    "BackPatchName": "",
                    "AirfoilPatchName": ""
                }
            }
        }
