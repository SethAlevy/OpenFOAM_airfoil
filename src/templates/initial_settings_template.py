from abc import ABC, abstractmethod


class Settings(ABC):

    @property
    @abstractmethod
    def airfoil_settings(self) -> dict:
        pass

    @property
    @abstractmethod
    def mesh_settings(self) -> dict:
        pass
