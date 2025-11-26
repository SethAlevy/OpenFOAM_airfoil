import numpy as np
from abc import ABC, abstractmethod


class Airfoil(ABC):
    """ 
    Base class for airfoil representation. Defines the essential properties that any
    airfoil class must implement.
    """

    @property
    @abstractmethod
    def upper_surface(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Get the upper surface coordinates eventually rotated by alpha.

        Returns:
            tuple[np.ndarray, np.ndarray]: Two arrays with x and y coordinates.
        """
        pass

    @property
    @abstractmethod
    def lower_surface(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Get the lower surface coordinates eventually rotated by alpha.

        Returns:
            tuple[np.ndarray, np.ndarray]: Two arrays with x and y coordinates.
        """
        pass

    @property
    @abstractmethod
    def mean_camber_line(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Get the mean camber line coordinates eventually rotated by alpha.

        Returns:
            tuple[np.ndarray, np.ndarray]: Two arrays with x and y coordinates.
        """
        pass

    @property
    @abstractmethod
    def thickness(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Get the thickness distribution coordinates eventually rotated by alpha.

        Returns:
            tuple[np.ndarray, np.ndarray]: Two arrays with x and y coordinates.
        """
        pass

    @property
    @abstractmethod
    def chord(self) -> float:
        """
        Get the chord length.

        Returns:
            float: chord length.
        """
        pass

    @property
    @abstractmethod
    def alpha(self) -> float:
        """
        Get the angle of attack.

        Returns:
            float: angle of attack in degrees.
        """
        pass
