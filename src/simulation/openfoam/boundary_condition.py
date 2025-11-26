from pathlib import Path
import numpy as np
from templates.initial_settings_template import Settings
import utils.utilities as ut
from utils.logger import SimpleLogger


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
            kinematic_pressure: bool = None,
            setup: Settings = None
    ) -> None:
        """
        Initialize boundary conditions based on provided parameters or setup settings
        if parameters are not provided.

        There are three main ways to define the boundary conditions:
        1. Directly provide the velocity.
        2. Provide the Mach number along with other parameters to calculate velocity.
        3. Provide the Reynolds number along with other parameters to calculate
        velocity.

        If there are more than one way to define the boundary conditions, the priority
        is same as above. The set of required parameters is not defined directly, as
        there are different combinations possible. Missing parameters will be calculated
        using general equations and models (like ISA for atmospheric conditions), but
        it is recommended to carefully define the desired initial conditions and avoid
        passing too many parameters at once.

        Args:
            velocity (np.ndarray): Velocity vector.
            mach_number (float): Mach number.
            reynolds_number (float): Reynolds number.
            density (float): Density value.
            temperature (float): Temperature value.
            altitude (float): Altitude value.
            nu (float): Kinematic viscosity.
            pressure (float): Pressure value.
            chord_length (float): Chord length of the airfoil.
            kinematic_pressure (bool): Flag indicating if pressure is kinematic.
            setup (Settings): Settings object containing simulation and airfoil settings.
        """

        self.bc_setup = setup.simulation_settings.get("BoundaryConditions", {})
        self.fluid_setup = setup.simulation_settings.get("Fluid", {})

        if kinematic_pressure is None:
            self._kinematic_pressure = self.bc_setup.get("KinematicPressure", None)
            if self._kinematic_pressure is None:
                self._kinematic_pressure = False
                SimpleLogger.warning(
                    "KinematicPressure flag not specified, assuming absolute pressure."
                )

        if velocity is None:
            self._velocity = self.bc_setup.get("Velocity", None)
        if mach_number is None:
            self._mach_number = self.bc_setup.get("MachNumber", None)
        if reynolds_number is None:
            self._reynolds_number = self.bc_setup.get("ReynoldsNumber", None)
        if density is None:
            self._density = self.fluid_setup.get("Density", None)
        if temperature is None:
            self._temperature = self.fluid_setup.get("Temperature", None)
        if altitude is None:
            self._altitude = self.fluid_setup.get("Altitude", None)
        if nu is None:
            self._nu = self.fluid_setup.get("KinematicViscosity", None)
        if pressure is None:
            self._pressure = self.bc_setup.get("Pressure", None)
        if chord_length is None:
            chord_length = setup.airfoil_settings.get("Chord", None)
            if chord_length is None:
                raise ValueError(
                    "Chord length must be provided either directly or in the setup.")
            else:
                self._chord = chord_length

        if altitude is not None and (
                temperature is None or pressure is None or density is None):
            SimpleLogger.warning(
                "Altitude provided, calculating temperature, pressure, and density"
                " using ISA if missing."
            )

            if pressure is None and self._kinematic_pressure:
                convert_to_kinematic = True

            temperature, pressure, density = ut.international_standard_atmosphere(
                altitude, temperature, pressure, density)

            SimpleLogger.log(
                f"Conditions: {altitude} m, temperature: {temperature} K,"
                f" absolute pressure: {pressure} Pa, ρ={density} kg/m³"
            )

            if convert_to_kinematic:
                pressure *= density
                SimpleLogger.warning(
                    "Converting calculated absolute pressure to kinematic pressure."
                )

            self._temperature = temperature
            self._pressure = pressure
            self._density = density

        if self.nu is None and self.temperature is not None and self.density is not None:
            SimpleLogger.warning(
                "Calculating kinematic viscosity based on temperature and density."
            )

            self._nu = ut.kinematic_viscosity_air(self.temperature, self.density)

            SimpleLogger.log(f"Kinematic viscosity: {self.nu} m²/s")

        if self.velocity is not None:
            self._as_velocity = True
            self.from_velocity(self.velocity, self.pressure)
        elif self.mach_number is not None and self.pressure is not None and self.temperature is not None:
            self._as_mach = True
            self.from_mach_number(self.mach_number, self.pressure, self.temperature)
        elif self.reynolds_number is not None and self.nu is not None and self.pressure is not None:
            self._as_reynolds = True
            self.from_reynolds_number(self.reynolds_number, self.nu, self.pressure)
        else:
            SimpleLogger.warning(
                "Insufficient parameters to determine boundary conditions. No"
                " parameters initialized automatically."
            )
        self.bc_details()

    def from_velocity(
            self,
            velocity: np.ndarray = None,
            pressure: float = None,
    ) -> None:
        """
        Set boundary conditions directly from provided velocity and pressure.

        Args:
            velocity (np.ndarray): Velocity vector.
            pressure (float): Pressure value.
        """
        SimpleLogger.log(
            "Velocity provided, setting boundary conditions from it (default variant)."
        )
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
        """
        Calculate velocity from Mach number if it is not provided directly and
        set boundary conditions.
        
        Args:
            mach_number (float): Mach number.
            pressure (float): Pressure value.
            temperature (float): Temperature value.
        """
        SimpleLogger.log(
            "Mach number provided, calculating velocity and setting boundary"
            " conditions."
        )
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
            reynolds_number: float = None,
            nu: float = None,
            pressure: float = None,
    ) -> None:
        """
        Calculate velocity from Reynolds number if it is not provided directly and set
        boundary conditions
        
        Args:
            reynolds_number (float): Reynolds number.
            nu (float): Kinematic viscosity.
            pressure (float): Pressure value.
        """
        SimpleLogger.log(
            "Reynolds number provided, calculating velocity and setting boundary"
            " conditions."
        )
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
        Generate the formatted string for the '0/U' boundary condition file with custom
        patch names and defined velocity compatible with OpenFOAM. To use custom
        values use set_velocity() method.

        Args:
            inlet_patch (str): Name of the inlet patch.
            outlet_patch (str): Name of the outlet patch.
            lower_wall_patch (str): Name of the lower wall patch.
            upper_wall_patch (str): Name of the upper wall patch.
            front_back_patch (str): Name of the front and back patch.

        Returns:
            str: Formatted string for the '0/U' boundary condition file.
        """
        velocity = self.velocity
        if isinstance(velocity, (list, tuple, np.ndarray)):
            velocity_str = " ".join(str(v) for v in velocity)
        else:
            velocity_str = f"{velocity} 0 0"

        if inlet_patch is None:
            inlet_patch = self.bc_setup.get("InletPatch", "inlet")
        if outlet_patch is None:
            outlet_patch = self.bc_setup.get("OutletPatch", "outlet")
        if lower_wall_patch is None:
            lower_wall_patch = self.bc_setup.get(
                "LowerWallPatch", "lowerWall")
        if upper_wall_patch is None:
            upper_wall_patch = self.bc_setup.get(
                "UpperWallPatch", "upperWall")
        if front_back_patch is None:
            front_back_patch = self.bc_setup.get(
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
    
    def to_absolute_pressure(self) -> None:
        """
        Convert kinematic pressure to absolute pressure.

        Args:
            pressure (float): Gauge pressure value.

        Returns:
            float: Absolute pressure value.
        """
        self._pressure /= self.density
        self._kinematic_pressure = False

    def to_kinematic_pressure(self) -> None:
        """
        Convert absolute pressure to kinematic pressure.

        Args:
            pressure (float): Absolute pressure value.

        Returns:
            float: Kinematic pressure value.
        """
        self._pressure *= self.density
        self._kinematic_pressure = True

    def pressure_bc(
            self,
            inlet_patch: str = None,
            outlet_patch: str = None,
            lower_wall_patch: str = None,
            upper_wall_patch: str = None,
            front_back_patch: str = None
    ) -> str:
        """
        Generate the formatted string for the '0/p' boundary condition file with custom
        patch names and defined pressure compatible with OpenFOAM. To use custom
        values use set_pressure() method.

        Args:
            inlet_patch (str): Name of the inlet patch.
            outlet_patch (str): Name of the outlet patch.
            lower_wall_patch (str): Name of the lower wall patch.
            upper_wall_patch (str): Name of the upper wall patch.
            front_back_patch (str): Name of the front and back patch.

        Returns:
            str: Formatted string for the '0/p' boundary condition file.
        """
        if inlet_patch is None:
            inlet_patch = self.bc_setup.get("InletPatch", "inlet")
        if outlet_patch is None:
            outlet_patch = self.bc_setup.get("OutletPatch", "outlet")
        if lower_wall_patch is None:
            lower_wall_patch = self.bc_setup.get("LowerWallPatch", "lowerWall")
        if upper_wall_patch is None:
            upper_wall_patch = self.bc_setup.get("UpperWallPatch", "upperWall")
        if front_back_patch is None:
            front_back_patch = self.bc_setup.get("FrontBackPatch", "frontAndBack")

        if self._kinematic_pressure:
            pressure_value = self.pressure
        else:
            SimpleLogger.warning(
                "Pressure boundary condition set as absolute pressure, converting"
                " to kinematic for OpenFOAM."
            )
            pressure_value = self.pressure / self.density

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
    
    def set_custom_velocity(self, velocity: np.ndarray) -> None:
        """
        Set a custom velocity vector.

        Args:
            velocity (np.ndarray): The velocity vector to set.
        """
        self._velocity = velocity

    def set_custom_pressure(self, pressure: float, kinematic: bool = False) -> None:
        """
        Set a custom pressure value.

        Args:
            pressure (float): The pressure value to set.
            kinematic (bool): Flag indicating if the pressure is kinematic.
        """
        self._kinematic_pressure = kinematic

        if kinematic:
            self._pressure = pressure * self.density
        else:
            self._pressure = pressure

    def write_bc(self, content: str, output_path: Path) -> None:
        """
        Write the velocity or pressure boundary condition content to a file.

        Args:
            content (str): The content to write to the file.
            output_path (Path): The path to the output file.
        """
        with open(output_path, 'w') as file:
            file.write(content)

    def bc_details(self) -> None:
        """
        Log the details of the current boundary conditions.
        """
        SimpleLogger.log("Boundary Conditions Details:")
        if self._as_velocity:
            SimpleLogger.log("Boundary conditions define directly from velocity")
        elif self._as_mach:
            SimpleLogger.log("Boundary conditions define from Mach number")
        elif self._as_reynolds:
            SimpleLogger.log("Boundary conditions define from Reynolds number")
        
        SimpleLogger.log(
            "Main parameters: \n"
            f"Velocity: {self.velocity} m/s\n"
        )
        if self._kinematic_pressure:
            SimpleLogger.log(f"Pressure: {self.pressure} (kinematic) m²/s²\n")
        else:
            SimpleLogger.log(f"Pressure: {self.pressure} Pa\n")

        SimpleLogger.log(
            "Additional parameters: \n"
            f"Mach Number: {self.mach_number}\n"
            f"Reynolds Number: {self.reynolds_number}\n"
            f"Altitude: {self.altitude} m\n"
        )
        SimpleLogger.log(
            "Fluid parameters: \n"
            f"Kinematic Viscosity: {self.nu} m²/s"
            f"Density: {self.density} kg/m³\n"
            f"Temperature: {self.temperature} K\n"
        )
        SimpleLogger.log(f"Chord Length: {self.chord} m")

    @property
    def velocity(self) -> np.ndarray:
        """
        Get the current velocity vector.

        Returns:
            np.ndarray: The current velocity vector.
        """
        return self._velocity

    @property
    def pressure(self) -> float:
        """
        Get the current pressure value.

        Returns:
            float: The current pressure value.
        """
        return self._pressure

    @property
    def mach_number(self) -> float:
        """
        Get the current Mach number.

        Returns:
            float: The current Mach number.
        """
        return self._mach_number

    @property
    def reynolds_number(self) -> float:
        """
        Get the current Reynolds number.

        Returns:
            float: The current Reynolds number.
        """
        return self._reynolds_number

    @property
    def density(self) -> float:
        """
        Get the current density value.

        Returns:
            float: The current density value.
        """
        return self._density

    @property
    def temperature(self) -> float:
        """
        Get the current temperature value.

        Returns:
            float: The current temperature value.
        """
        return self._temperature

    @property
    def altitude(self) -> float:
        """
        Get the current altitude value.

        Returns:
            float: The current altitude value.
        """
        return self._altitude

    @property
    def nu(self) -> float:
        """
        Get the current kinematic viscosity value.

        Returns:
            float: The current kinematic viscosity value.
        """
        return self._nu

    @property
    def chord(self) -> float:
        """
        Get the chord length.

        Returns:
            float: The chord length.
        """
        return self._chord
