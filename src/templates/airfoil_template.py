import numpy as np
from abc import ABC, abstractmethod


class Airfoil(ABC):

    @property
    @abstractmethod
    def upper_surface(self) -> tuple[np.ndarray, np.ndarray]:
        pass

    @property
    @abstractmethod
    def lower_surface(self) -> tuple[np.ndarray, np.ndarray]:
        pass

    @property
    @abstractmethod
    def mean_camber_line(self) -> tuple[np.ndarray, np.ndarray]:
        pass

    @property
    @abstractmethod
    def thickness(self) -> tuple[np.ndarray, np.ndarray]:
        pass

    @property
    @abstractmethod
    def chord(self) -> float:
        pass

    @property
    @abstractmethod
    def alpha(self) -> float:
        pass
