import numpy as np
from abc import ABC, abstractmethod


class BoundaryCondition(ABC):
    """ 
    Base class for boundary condition representation. Defines the essential properties that any
    boundary condition class must implement.
    """

    @property
    @abstractmethod
    def velocity(self) -> np.ndarray:
        """
        Get the current velocity vector.

        Returns:
            np.ndarray: The current velocity vector.
        """
        pass
    
    @property
    @abstractmethod
    def pressure(self) -> float:
        """
        Get the current pressure value.

        Returns:
            float: The current pressure value.
        """
        pass

    @property
    @abstractmethod
    def mach_number(self) -> float:
        """
        Get the current Mach number.

        Returns:
            float: The current Mach number.
        """
        pass

    @property
    @abstractmethod
    def reynolds_number(self) -> float:
        """
        Get the current Reynolds number.

        Returns:
            float: The current Reynolds number.
        """
        pass

    @property
    @abstractmethod
    def density(self) -> float:
        """
        Get the current density value.

        Returns:
            float: The current density value.
        """
        pass

    @property
    @abstractmethod
    def temperature(self) -> float:
        """
        Get the current temperature value.

        Returns:
            float: The current temperature value.
        """
        pass

    @property
    @abstractmethod
    def altitude(self) -> float:
        """
        Get the current altitude value.

        Returns:
            float: The current altitude value.
        """
        pass

    @property
    @abstractmethod
    def nu(self) -> float:
        """
        Get the current kinematic viscosity value.

        Returns:
            float: The current kinematic viscosity value.
        """
        pass

    @property
    @abstractmethod
    def chord(self) -> float:
        """
        Get the chord length.

        Returns:
            float: The chord length.
        """
        pass
