import numpy as np
from abc import ABC, abstractmethod
from pathlib import Path


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

    @abstractmethod
    def to_stl(self, output_path: Path, thickness: float, dimension: int) -> None:
        """
        Export airfoil geometry to STL format. It can be 2D or 3D by selecting the 
        dimension parameter (in general depends on the mesh requirements).

        Args:
            output_path (Path): Path to save the STL file.
            thickness (float): Thickness for 3D extrusion (only used if dimension=3).
            dimension (int): 2 for 2D STL, 3 for extruded 3D STL.
        """
        pass

    @abstractmethod
    def plot(
        self,
        title: str = "Airfoil",
        save_path: Path = None,
        show: bool = True
    ) -> None:
        """
        Plot the airfoil geometry.

        Args:
            title (str): Title of the plot.
            save_path (Path): Path to save the plot image. If None, the plot is shown.
            show (bool): Whether to display the plot.
        """
        pass

    @abstractmethod
    def airfoil_details(self) -> None:
        """
        Simple log with the airfoil details.
        """
        pass
