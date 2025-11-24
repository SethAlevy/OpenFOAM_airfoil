from pathlib import Path
import numpy as np
from templates.initial_settings_template import Settings
import utils.utilities as ut


class BoundaryConditions:
    def __init__(
            self,
            velocity: np.ndarray = None,
            mach_number: float = None,
            reynolds_number: float = None,
            density: float = None,
            temperature: float = None,
            altitude: float = None,
            nu: float = None,
            pressure: float = None,
            chord_length: float = None,
            setup: Settings = None
    ) -> None:

        self.boundary_conditions_setup = setup.simulation_settings.get("BoundaryConditions", {})
        self.fluid_setup = setup.simulation_settings.get("Fluid", {})

        if velocity is None:
            self._velocity = self.boundary_conditions_setup.get("Velocity", None)
        if mach_number is None:
            self._mach_number = self.boundary_conditions_setup.get("MachNumber", None)
        if reynolds_number is None:
            self._reynolds_number = self.boundary_conditions_setup.get(
                "ReynoldsNumber", None)
        if density is None:
            self._density = self.fluid_setup.get("Density", None)
        if temperature is None:
            self._temperature = self.fluid_setup.get("Temperature", None)
        if altitude is None:
            self._altitude = self.fluid_setup.get("Altitude", None)
        if nu is None:
            self._nu = self.fluid_setup.get("KinematicViscosity", None)
        if pressure is None:
            self._pressure = self.boundary_conditions_setup.get("Pressure", None)
        if chord_length is None:
            chord_length = setup.airfoil_settings.get("Chord", None)
            if chord_length is None:
                raise ValueError(
                    "Chord length must be provided either directly or in the setup.")
            else:
                self._chord = chord_length

        if altitude is not None:
            print(
                "Altitude provided, calculating temperature, pressure, and density using ISA if missing.")
            temperature, pressure, density = ut.international_standard_atmosphere(
                altitude, temperature, pressure, density)
            print(
                f"Conditions: {altitude} m, T={temperature} K, p={pressure} Pa, ρ={density} kg/m³")

            self._temperature = temperature
            self._pressure = pressure
            self._density = density

        if self.nu is None and self.temperature is not None and self.density is not None:
            print("Calculating kinematic viscosity based on temperature and density.")
            self._nu = ut.kinematic_viscosity_air(self.temperature, self.density)
            print(f"Kinematic viscosity: {self.nu} m²/s")

        if self.velocity is not None:
            self.from_velocity(self.velocity, self.pressure)
        elif self.mach_number is not None and self.pressure is not None and self.temperature is not None:
            self.from_mach_number(self.mach_number, self.pressure, self.temperature)
        elif self.reynolds_number is not None and self.nu is not None and self.pressure is not None:
            self.from_reynolds_number(self.reynolds_number, self.nu, self.pressure)
        else:
            print(
                "Insufficient parameters to determine boundary conditions. No parameters" \
                "initialized automatically.")

    def from_velocity(
            self,
            velocity: np.ndarray = None,
            pressure: float = None,
    ) -> None:
        print("Velocity provided, setting boundary conditions from it (default variant).")
        if velocity is not None:
            self._velocity = velocity
        if pressure is not None:
            self._pressure = pressure

    def from_mach_number(
            self,
            mach_number: float = None,
            pressure: float = None,
            temperature: float = None,
    ) -> None:
        print("Mach number provided, calculating velocity and setting boundary conditions.")
        if mach_number is not None:
            self._mach_number = mach_number
        if pressure is not None:
            self._pressure = pressure
        if temperature is not None:
            self._temperature = temperature

        speed_of_sound = np.sqrt(1.4 * 287.05 * self.temperature)
        self._velocity = mach_number * speed_of_sound

    def from_reynolds_number(
            self,
            reynolds_number: float,
            nu: float,
            pressure: float,
    ) -> None:
        print("Reynolds number provided, calculating velocity and setting boundary conditions.")
        if reynolds_number is not None:
            self._reynolds_number = reynolds_number
        if nu is not None:
            self._nu = nu
        if pressure is not None:
            self._pressure = pressure

        self._velocity = self.reynolds_number * self.nu / self.chord

    def velocity_bc(
            self,
            inlet_patch: str = None,
            outlet_patch: str = None,
            lower_wall_patch: str = None,
            upper_wall_patch: str = None,
            front_back_patch: str = None
    ) -> str:
        """
        Generate the formatted string for the '0/U' boundary condition file with custom patch names.
        """
        velocity = self.velocity
        if isinstance(velocity, (list, tuple, np.ndarray)):
            velocity_str = " ".join(str(v) for v in velocity)
        else:
            velocity_str = f"{velocity} 0 0"

        if inlet_patch is None:
            inlet_patch = self.boundary_conditions_setup.get("InletPatch", "inlet")
        if outlet_patch is None:
            outlet_patch = self.boundary_conditions_setup.get("OutletPatch", "outlet")
        if lower_wall_patch is None:
            lower_wall_patch = self.boundary_conditions_setup.get(
                "LowerWallPatch", "lowerWall")
        if upper_wall_patch is None:
            upper_wall_patch = self.boundary_conditions_setup.get(
                "UpperWallPatch", "upperWall")
        if front_back_patch is None:
            front_back_patch = self.boundary_conditions_setup.get(
                "FrontBackPatch", "frontAndBack")

        return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform ({velocity_str});

boundaryField
{{
    {inlet_patch}
    {{
        type            fixedValue;
        value           uniform ({velocity_str});
    }}
    {outlet_patch}
    {{
        type            zeroGradient;
    }}
    {lower_wall_patch}
    {{
        type            noSlip;
    }}
    {upper_wall_patch}
    {{
        type            noSlip;
    }}
    {front_back_patch}
    {{
        type            empty;
    }}
}}
"""

    def pressure_bc(
            self,
            inlet_patch: str = None,
            outlet_patch: str = None,
            lower_wall_patch: str = None,
            upper_wall_patch: str = None,
            front_back_patch: str = None
    ) -> str:
        """
        Generate the formatted string for the '0/p' boundary condition file with custom patch names.
        """
        if inlet_patch is None:
            inlet_patch = self.boundary_conditions_setup.get("InletPatch", "inlet")
        if outlet_patch is None:
            outlet_patch = self.boundary_conditions_setup.get("OutletPatch", "outlet")
        if lower_wall_patch is None:
            lower_wall_patch = self.boundary_conditions_setup.get(
                "LowerWallPatch", "lowerWall")
        if upper_wall_patch is None:
            upper_wall_patch = self.boundary_conditions_setup.get(
                "UpperWallPatch", "upperWall")
        if front_back_patch is None:
            front_back_patch = self.boundary_conditions_setup.get(
                "FrontBackPatch", "frontAndBack")

        pressure_value = self.pressure

        return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform {pressure_value};

boundaryField
{{
    {inlet_patch}
    {{
        type            zeroGradient;
    }}
    {outlet_patch}
    {{
        type            fixedValue;
        value           uniform {pressure_value};
    }}
    {lower_wall_patch}
    {{
        type            zeroGradient;
    }}
    {upper_wall_patch}
    {{
        type            zeroGradient;
    }}
    {front_back_patch}
    {{
        type            empty;
    }}
}}
"""

    def write_bc(self, content: str, output_path: Path) -> None:
        """
        Write the velocity or pressure boundary condition content to a file.

        Args:
            content (str): The content to write to the file.
            output_path (Path): The path to the output file.
        """
        with open(output_path, 'w') as file:
            file.write(content)

    @property
    def velocity(self) -> np.ndarray:
        return self._velocity

    @property
    def pressure(self) -> float:
        return self._pressure

    @property
    def mach_number(self) -> float:
        return self._mach_number

    @property
    def reynolds_number(self) -> float:
        return self._reynolds_number

    @property
    def density(self) -> float:
        return self._density

    @property
    def temperature(self) -> float:
        return self._temperature

    @property
    def altitude(self) -> float:
        return self._altitude

    @property
    def nu(self) -> float:
        return self._nu

    @property
    def chord(self) -> float:
        return self._chord
