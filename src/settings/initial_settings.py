import json


class InitialSettingsReader():
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.settings = self.read_settings()

    def read_settings(self) -> dict:
        settings = {}
        with open(self.filepath, 'r') as file:
            settings = json.load(file)
        return settings
    
    @property
    def airfoil_settings(self) -> dict:
        return self.settings.get("Airfoil", {})

    @property
    def mesh_settings(self) -> dict:
        return self.settings.get("Mesh", {})
