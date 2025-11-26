from abc import ABC, abstractmethod


class Settings(ABC):

    @property
    @abstractmethod
    def airfoil_settings(self) -> dict:
        """
        Contains parameters related to airfoil generation, like designation, chord.

        Returns:
            dict: Airfoil parameters.
        """
        pass

    @property
    @abstractmethod
    def mesh_settings(self) -> dict:
        """
        Contains parameters related to mesh generation in openFOAM. It allows to handle
        the creation of blockMesh and further refinement using snappyHexMesh or cfMesh.

        Returns:
            dict: Mesh parameters.
        """
        pass

    @property
    @abstractmethod
    def simulation_settings(self) -> dict:
        """
        Contains parameters to handle the physical properties of the simulation 
        like turbulence model, fluid properties and openFOAM files like controlDict,
        fvSolution, fvSchemes and boundary conditions.

        Returns:
            dict: Simulation parameters.
        """
        pass

    @staticmethod
    @abstractmethod
    def template_settings() -> dict:
        """
        A full settings template to allow validation checks on the settings json.

        Returns:
            dict: Template settings dictionary.
        """
        pass
