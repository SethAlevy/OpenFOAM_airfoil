from pathlib import Path
import numpy as np
from templates.initial_settings_template import Settings
import utils.utilities as ut
from utils.logger import SimpleLogger
import csv


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
            kinematic_pressure: bool = None,
            turbulence_model: str = None,
            turbulence_intensity: float = None,
            turbulence_length_scale: float = None,
            chord_length: float = None,
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
            kinematic_pressure (bool): Flag indicating if pressure is kinematic.
            turbulence_model (str): Selected turbulence model. If not provided, it will
            be taken from the setup settings or assumed as "kOmegaSST".
            turbulence_intensity (float): Turbulence intensity.
            turbulence_length_scale (float): Turbulence length scale.
            chord_length (float): Chord length of the airfoil.
            setup (Settings): Settings object containing simulation and airfoil
                settings.
        """

        self.bc_setup = setup.simulation_settings.get("BoundaryConditions", {})
        self.fluid_setup = setup.simulation_settings.get("Fluid", {})
        self.turbulence_setup = setup.simulation_settings.get(
            "TurbulenceProperties", {})

        if kinematic_pressure is None:
            self._kinematic_pressure = self.bc_setup.get("KinematicPressure", None)
            if self._kinematic_pressure is None:
                self._kinematic_pressure = False
                SimpleLogger.warning(
                    "KinematicPressure flag not specified, assuming absolute pressure."
                )

        if velocity is None:
            self._velocity = self.bc_setup.get("Velocity", None)
        else:
            self._velocity = velocity
        if mach_number is None:
            self._mach_number = self.bc_setup.get("MachNumber", None)
        else:
            self._mach_number = mach_number
        if reynolds_number is None:
            self._reynolds_number = self.bc_setup.get("ReynoldsNumber", None)
        else:
            self._reynolds_number = reynolds_number
        if density is None:
            self._density = self.fluid_setup.get("Density", None)
        else:
            self._density = density
        if temperature is None:
            self._temperature = self.fluid_setup.get("Temperature", None)
        else:
            self._temperature = temperature
        if altitude is None:
            self._altitude = self.fluid_setup.get("Altitude", None)
        else:
            self._altitude = altitude
        if nu is None:
            self._nu = self.fluid_setup.get("KinematicViscosity", None)
        else:
            self._nu = nu
        if pressure is None:
            self._pressure = self.bc_setup.get("Pressure", None)
        else:
            self._pressure = pressure

        if turbulence_model is None:
            turbulence_model = self.turbulence_setup.get("Model", None)
            if turbulence_model is None:
                SimpleLogger.warning(
                    "Turbulence model no specified, assuming default kOmegaSST."
                )
                turbulence_model = "kOmegaSST"
                turbulence_intensity = 0.01
                turbulence_length_scale = 0.08
        self._turbulence_model = turbulence_model
        if turbulence_intensity is None:
            self._turbulence_intensity = self.turbulence_setup.get(
                "TurbulenceIntensity", None)
        else:
            self.turbulence_setup["TurbulenceIntensity"] = turbulence_intensity
        if turbulence_length_scale is None:
            self._turbulence_length_scale = self.turbulence_setup.get(
                "TurbulenceLengthScale", None)
        else:
            self.turbulence_setup["TurbulenceLengthScale"] = turbulence_length_scale

        if chord_length is None:
            chord_length = setup.airfoil_settings.get("Chord", None)
            if chord_length is None:
                raise ValueError(
                    "Chord length must be provided either directly or in the setup.")
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

        if (
            self.nu is None and self.temperature is not None
            and self.density is not None
        ):
            SimpleLogger.warning(
                "Calculating kinematic viscosity based on temperature and density."
            )

            self._nu = ut.kinematic_viscosity_air(self.temperature, self.density)

            SimpleLogger.log(f"Kinematic viscosity: {self.nu} m²/s")

        self._as_velocity = False
        self._as_mach = False
        self._as_reynolds = False

        if self.velocity is not None:
            self._as_velocity = True
            self.from_velocity(self.velocity, self.pressure)
        elif (
            self.mach_number is not None and self.pressure is not None
            and self.temperature is not None
        ):
            self._as_mach = True
            self.from_mach_number(self.mach_number, self.pressure, self.temperature)
        elif (
            self.reynolds_number is not None and self.nu is not None
            and self.pressure is not None
        ):
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

        reynolds_number = (self.velocity * self.chord) / self.nu
        self._reynolds_number = reynolds_number

        if self.temperature is not None:
            speed_of_sound = np.sqrt(1.4 * 287.05 * self.temperature)
            mach_number = np.linalg.norm(self.velocity) / speed_of_sound
            self._mach_number = mach_number
        else:
            self._mach_number = None

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

        reynolds_number = (self.velocity * self.chord) / self.nu
        self._reynolds_number = reynolds_number

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

        if self.temperature is not None:
            speed_of_sound = np.sqrt(1.4 * 287.05 * self.temperature)
            mach_number = np.linalg.norm(self.velocity) / speed_of_sound
            self._mach_number = mach_number
        else:
            self._mach_number = None

    def velocity_bc(
            self,
            inlet_patch: str = None,
            outlet_patch: str = None,
            lower_wall_patch: str = None,
            upper_wall_patch: str = None,
            front_patch: str = None,
            back_patch: str = None,
            airfoil_patch: str = None
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
            front_patch (str): Name of the front patch.
            back_patch (str): Name of the back patch.
            airfoil_patch (str): Name of the airfoil patch.

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
            lower_wall_patch = self.bc_setup.get("LowerWallPatch", "lowerWall")
        if upper_wall_patch is None:
            upper_wall_patch = self.bc_setup.get("UpperWallPatch", "upperWall")
        if front_patch is None:
            front_patch = self.bc_setup.get("FrontPatch", "front")
        if back_patch is None:
            back_patch = self.bc_setup.get("BackPatch", "back")
        if airfoil_patch is None:
            airfoil_patch = self.bc_setup.get("AirfoilPatch", "airfoil")

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
        type            zeroGradient;  // or noSlip if this is a wall
    }}
    {upper_wall_patch}
    {{
        type            zeroGradient;  // or noSlip if this is a wall
    }}
    {front_patch}
    {{
        type            symmetryPlane;
    }}
    {back_patch}
    {{
        type            symmetryPlane;
    }}
    {airfoil_patch}
    {{
        type            noSlip;
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
            front_patch: str = None,
            back_patch: str = None,
            airfoil_patch: str = None
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
            front_patch (str): Name of the front patch.
            back_patch (str): Name of the back patch.
            airfoil_patch (str): Name of the airfoil patch.

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
        if front_patch is None:
            front_patch = self.bc_setup.get("FrontPatch", "front")
        if back_patch is None:
            back_patch = self.bc_setup.get("BackPatch", "back")
        if airfoil_patch is None:
            airfoil_patch = self.bc_setup.get("AirfoilPatch", "airfoil")

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
    {front_patch}
    {{
        type            symmetryPlane;
    }}
    {back_patch}
    {{
        type            symmetryPlane;
    }}
    {airfoil_patch}
    {{
        type            zeroGradient;
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

    def turbulence_bc(self) -> dict[str, str]:
        """
        Generate the formatted string for the turbulence boundary conditions files for
        the current turbulence model. Supported models: kOmegaSST, kEpsilon and
        SpalartAllmaras.

        Returns:
            dict[str, str]: Formatted strings for the turbulence boundary condition
                files.
        """

        model = self._turbulence_model.lower()
        if model == "komegasst":
            return self.kOmegaSST_bc()
        elif model == "kepsilon":
            return self.kEpsilon_bc()
        elif model == "spalartallmaras":
            return self.SpalartAllmaras_bc()
        else:
            raise ValueError(f"Unsupported turbulence model: {model}")

    def kOmegaSST_bc(self,
                     inlet_patch: str = None,
                     outlet_patch: str = None,
                     lower_wall_patch: str = None,
                     upper_wall_patch: str = None,
                     front_patch: str = None,
                     back_patch: str = None,
                     airfoil_patch: str = None
                     ) -> dict[str, str]:
        """
        Generate boundary conditions for k-omega SST model (0/k, 0/omega, 0/nut).
        """
        if inlet_patch is None:
            inlet_patch = self.bc_setup.get("InletPatch", "inlet")
        if outlet_patch is None:
            outlet_patch = self.bc_setup.get("OutletPatch", "outlet")
        if lower_wall_patch is None:
            lower_wall_patch = self.bc_setup.get("LowerWallPatch", "lowerWall")
        if upper_wall_patch is None:
            upper_wall_patch = self.bc_setup.get("UpperWallPatch", "upperWall")
        if front_patch is None:
            front_patch = self.bc_setup.get("FrontPatch", "front")
        if back_patch is None:
            back_patch = self.bc_setup.get("BackPatch", "back")
        if airfoil_patch is None:
            airfoil_patch = self.bc_setup.get("AirfoilPatch", "airfoil")

        # 0/k
        k_bc = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform {self.k};

boundaryField
{{
    {inlet_patch}
    {{
        type            fixedValue;
        value           uniform {self.k};
    }}
    {outlet_patch}
    {{
        type            zeroGradient;
    }}
    {lower_wall_patch}
    {{
        type            zeroGradient;
    }}
    {upper_wall_patch}
    {{
        type            zeroGradient;
    }}
    {front_patch}
    {{
        type            symmetryPlane;
    }}
    {back_patch}
    {{
        type            symmetryPlane;
    }}
    {airfoil_patch}
    {{
        type            kqRWallFunction;
        value           uniform 0;
    }}
}}
"""

        # 0/omega
        omega_bc = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      omega;
}}

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform {self.omega};

boundaryField
{{
    {inlet_patch}
    {{
        type            fixedValue;
        value           uniform {self.omega};
    }}
    {outlet_patch}
    {{
        type            zeroGradient;
    }}
    {lower_wall_patch}
    {{
        type            zeroGradient;
    }}
    {upper_wall_patch}
    {{
        type            zeroGradient;
    }}
    {front_patch}
    {{
        type            symmetryPlane;
    }}
    {back_patch}
    {{
        type            symmetryPlane;
    }}
    {airfoil_patch}
    {{
        type            omegaWallFunction;
        value           uniform 0;
    }}
}}
"""

        # 0/nut
        nut_bc = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}}

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
    {inlet_patch}
    {{
        type            zeroGradient;
    }}
    {outlet_patch}
    {{
        type            zeroGradient;
    }}
    {lower_wall_patch}
    {{
        type            zeroGradient;
    }}
    {upper_wall_patch}
    {{
        type            zeroGradient;
    }}
    {front_patch}
    {{
        type            symmetryPlane;
    }}
    {back_patch}
    {{
        type            symmetryPlane;
    }}
    {airfoil_patch}
    {{
        type            nutUSpaldingWallFunction;
        value           uniform 0;
    }}
}}
"""

        return {"k": k_bc, "omega": omega_bc, "nut": nut_bc}

    def kEpsilon_bc(self,
                    inlet_patch: str = None,
                    outlet_patch: str = None,
                    lower_wall_patch: str = None,
                    upper_wall_patch: str = None,
                    front_patch: str = None,
                    back_patch: str = None,
                    airfoil_patch: str = None
                    ) -> dict[str, str]:
        """
        Generate boundary conditions for k-epsilon model (0/k, 0/epsilon, 0/nut).
        """
        if inlet_patch is None:
            inlet_patch = self.bc_setup.get("InletPatch", "inlet")
        if outlet_patch is None:
            outlet_patch = self.bc_setup.get("OutletPatch", "outlet")
        if lower_wall_patch is None:
            lower_wall_patch = self.bc_setup.get("LowerWallPatch", "lowerWall")
        if upper_wall_patch is None:
            upper_wall_patch = self.bc_setup.get("UpperWallPatch", "upperWall")
        if front_patch is None:
            front_patch = self.bc_setup.get("FrontPatch", "front")
        if back_patch is None:
            back_patch = self.bc_setup.get("BackPatch", "back")
        if airfoil_patch is None:
            airfoil_patch = self.bc_setup.get("AirfoilPatch", "airfoil")

        # 0/k
        k_bc = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform {self.k};

boundaryField
{{
    {inlet_patch}
    {{
        type            fixedValue;
        value           uniform {self.k};
    }}
    {outlet_patch}
    {{
        type            zeroGradient;
    }}
    {lower_wall_patch}
    {{
        type            zeroGradient;
    }}
    {upper_wall_patch}
    {{
        type            zeroGradient;
    }}
    {front_patch}
    {{
        type            symmetryPlane;
    }}
    {back_patch}
    {{
        type            symmetryPlane;
    }}
    {airfoil_patch}
    {{
        type            fixedValue;
        value           uniform 0;
    }}
}}
"""

        # 0/epsilon
        epsilon_bc = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      epsilon;
}}

dimensions      [0 2 -3 0 0 0 0];

internalField   uniform {self.epsilon};

boundaryField
{{
    {inlet_patch}
    {{
        type            fixedValue;
        value           uniform {self.epsilon};
    }}
    {outlet_patch}
    {{
        type            zeroGradient;
    }}
    {lower_wall_patch}
    {{
        type            zeroGradient;
    }}
    {upper_wall_patch}
    {{
        type            zeroGradient;
    }}
    {front_patch}
    {{
        type            symmetryPlane;
        value           uniform {self.epsilon};
    }}
    {back_patch}
    {{
        type            symmetryPlane;
        value           uniform {self.epsilon};
    }}
    {airfoil_patch}
    {{
        type            fixedValue;
        value           uniform 0;
    }}
}}
"""

        # 0/nut
        nut_bc = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}}

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
    {inlet_patch}
    {{
        type            zeroGradient;
    }}
    {outlet_patch}
    {{
        type            zeroGradient;
    }}
    {lower_wall_patch}
    {{
        type            zeroGradient;
    }}
    {upper_wall_patch}
    {{
        type            zeroGradient;
    }}
    {front_patch}
    {{
        type            symmetryPlane;
    }}
    {back_patch}
    {{
        type            symmetryPlane;
    }}
    {airfoil_patch}
    {{
        type            nutUSpaldingWallFunction;
        value           uniform 0;
    }}
}}
"""

        return {"k": k_bc, "epsilon": epsilon_bc, "nut": nut_bc}

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
            f"                  Velocity: {self.velocity} m/s\n"
        )
        if self._kinematic_pressure:
            SimpleLogger.log(f"Pressure: {self.pressure} (kinematic) m²/s²\n")
        else:
            SimpleLogger.log(f"Pressure: {self.pressure} Pa\n")

        SimpleLogger.log(
            "Additional parameters: \n"
            f"                  Mach Number: {self.mach_number}\n"
            f"                  Reynolds Number: {self.reynolds_number}\n"
            f"                  Altitude: {self.altitude} m\n"
        )
        SimpleLogger.log(
            "Fluid parameters: \n"
            f"                  Kinematic Viscosity: {self.nu} m²/s \n"
            f"                  Density: {self.density} kg/m³\n"
            f"                  Temperature: {self.temperature} K\n"
        )
        SimpleLogger.log(f"Chord Length: {self.chord} m")

    def export_bc_to_csv(self, output_path: Path) -> None:
        """
        Export boundary conditions details to a CSV file.

        Args:
            output_path (Path): The path to the output CSV file.
        """
        # Determine definition method
        if self._as_velocity:
            definition_method = "Velocity"
        elif self._as_mach:
            definition_method = "Mach Number"
        elif self._as_reynolds:
            definition_method = "Reynolds Number"
        else:
            definition_method = "Unknown"

        # Calculate velocity magnitude
        if isinstance(self.velocity, np.ndarray):
            velocity_magnitude = np.linalg.norm(self.velocity)
        else:
            velocity_magnitude = self.velocity

        # Prepare pressure info
        pressure_unit = "m²/s²" if self._kinematic_pressure else "Pa"
        pressure_value = self.pressure if self.pressure is not None else "N/A"

        # Determine turbulence parameter
        if self.turbulence_model.lower() == "komegasst":
            turb_param = "omega"
            turb_unit = "1/s"
            turb_value = self.omega
        elif self.turbulence_model.lower() == "kepsilon":
            turb_param = "epsilon"
            turb_unit = "m²/s³"
            turb_value = self.epsilon
        else:
            turb_param = "N/A"
            turb_unit = "-"
            turb_value = "N/A"

        # Create data rows
        data = [
            ["Parameter", "Unit", "Value"],
            ["Definition Method", "-", definition_method],
            ["Velocity", "m/s", velocity_magnitude],
            ["Pressure", pressure_unit, pressure_value],
            ["Mach Number", "-",
             self.mach_number if self.mach_number is not None else "N/A"],
            ["Reynolds Number", "-",
             self.reynolds_number if self.reynolds_number is not None else "N/A"],
            ["Altitude", "m", self.altitude if self.altitude is not None else "N/A"],
            ["Kinematic Viscosity", "m2/s", self.nu if self.nu is not None else "N/A"],
            ["Density", "kg/m3", self.density if self.density is not None else "N/A"],
            ["Temperature", "K",
             self.temperature if self.temperature is not None else "N/A"],
            ["Chord Length", "m", self.chord],
            ["Turbulence Model", "-", self.turbulence_model],
            ["Turbulence Intensity", "-", self.turbulence_intensity],
            ["Turbulence Length Scale", "m", self.turbulence_length_scale],
            ["k", "m2/s2", self.k],
            [turb_param, turb_unit, turb_value]
        ]

        # Write to CSV
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(data)

        SimpleLogger.log(f"Boundary conditions exported to: {output_path}")

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

    @property
    def turbulence_model(self) -> str:
        """
        Get the turbulence model.

        Returns:
            str: The turbulence model.
        """
        return self._turbulence_model

    @property
    def k(self) -> float:
        """
        Get the turbulence kinetic energy.

        Returns:
            float: The turbulence kinetic energy.
        """
        U_inf = np.linalg.norm(self.velocity)

        k = 1.5 * (U_inf * self.turbulence_intensity) ** 2
        return k

    @property
    def omega(self) -> float:
        """
        Get the specific dissipation rate for k-omega SST model.

        Returns:
            float: The specific dissipation rate.
        """
        U_inf = np.linalg.norm(self.velocity)
        C_mu = 0.09

        k = 1.5 * (U_inf * self.turbulence_intensity) ** 2
        omega = np.sqrt(k) / (C_mu ** 0.25 * self.turbulence_length_scale)
        return omega

    @property
    def epsilon(self) -> float:
        """
        Get the dissipation rate for k-epsilon model.

        Returns:
            float: The dissipation rate.
        """
        U_inf = np.linalg.norm(self.velocity)
        C_mu = 0.09

        k = 1.5 * (U_inf * self.turbulence_intensity) ** 2
        epsilon = (C_mu ** 0.75) * (k ** 1.5) / self.turbulence_length_scale
        return epsilon

    @property
    def turbulence_intensity(self) -> float:
        """
        Get the turbulence intensity.

        Returns:
            float: The turbulence intensity.
        """
        return self._turbulence_intensity

    @property
    def turbulence_length_scale(self) -> float:
        """
        Get the turbulence length scale.

        Returns:
            float: The turbulence length scale.
        """
        return self._turbulence_length_scale
