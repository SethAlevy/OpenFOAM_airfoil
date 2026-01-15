import csv
import numpy as np
import utils.physics as ph
import utils.utilities as ut
from pathlib import Path
from utils.logger import SimpleLogger
from templates.initial_settings_template import Settings
from templates.openfoam_template_files.boundary_files import (
    U_bc,
    p_bc,
    k_bc,
    omega_bc,
    epsilon_bc,
    nut_bc,
    gammaInt_bc,
    retheta_bc
)


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
        Initialize boundary conditions based on provided parameters or setup
        settings if parameters are not provided.

        There are three main ways to define the boundary conditions:
        1. Directly provide the velocity.
        2. Provide the Mach number along with other parameters to calculate
           velocity.
        3. Provide the Reynolds number along with other parameters to calculate
           velocity.

        If there are more than one way to define the boundary conditions, the
        priority is same as above. The set of required parameters is not defined
        directly, as there are different combinations possible. Missing
        parameters will be calculated using general equations and models (like
        ISA for atmospheric conditions), but it is recommended to carefully
        define the desired initial conditions and avoid passing too many
        parameters at once.

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
            turbulence_model (str): Selected turbulence model. If not provided,
                it will be taken from the setup settings or assumed as
                "kOmegaSST".
            turbulence_intensity (float): Turbulence intensity.
            turbulence_length_scale (float): Turbulence length scale fraction. Needs
                multiplication by chord.
            chord_length (float): Chord length of the airfoil.
            setup (Settings): Settings object containing simulation and airfoil
                settings.
        """

        self.bc_setup = setup.simulation_settings.get("BoundaryConditions", {})
        self.fluid_setup = setup.simulation_settings.get("Fluid", {})
        self.turbulence_setup = setup.simulation_settings.get(
            "TurbulenceProperties", {})

        # Resolve kinematic_pressure flag
        self._kinematic_pressure = ut.resolve_value(
            kinematic_pressure, self.bc_setup, "KinematicPressure", False
        )
        if kinematic_pressure is None and "KinematicPressure" not in self.bc_setup:
            SimpleLogger.warning(
                "KinematicPressure flag not specified, assuming absolute "
                "pressure."
            )

        self._velocity = ut.resolve_value(velocity, self.bc_setup, "Velocity", None)
        self._mach_number = ut.resolve_value(
            mach_number, self.bc_setup, "MachNumber", None
        )
        self._reynolds_number = ut.resolve_value(
            reynolds_number, self.bc_setup, "ReynoldsNumber", None
        )
        self._density = ut.resolve_value(density, self.fluid_setup, "Density", None)
        self._temperature = ut.resolve_value(
            temperature, self.fluid_setup, "Temperature", None
        )
        self._altitude = ut.resolve_value(altitude, self.fluid_setup, "Altitude", None)
        self._nu = ut.resolve_value(nu, self.fluid_setup, "KinematicViscosity", None)
        self._pressure = ut.resolve_value(pressure, self.bc_setup, "Pressure", None)

        # Resolve turbulence model
        self._turbulence_model = ut.resolve_value(
            turbulence_model, self.turbulence_setup, "Model", None
        )
        if self._turbulence_model is None:
            SimpleLogger.warning(
                "Turbulence model not specified, assuming default "
                "kOmegaSST."
            )
            self._turbulence_model = "kOmegaSST"
            turbulence_intensity = 0.0003
            turbulence_length_scale = 0.015

        # Resolve turbulence parameters
        self._turbulence_intensity = ut.resolve_value(
            turbulence_intensity, self.turbulence_setup, "TurbulenceIntensity", None
        )
        if turbulence_intensity is not None:
            self.turbulence_setup["TurbulenceIntensity"] = turbulence_intensity

        self._turbulence_length_scale = ut.resolve_value(
            turbulence_length_scale,
            self.turbulence_setup,
            "TurbulenceLengthScaleChord",
            None
        )
        if turbulence_length_scale is not None:
            self.turbulence_setup["TurbulenceLengthScaleChord"] = (
                turbulence_length_scale
            )

        self._chord = ut.resolve_value(
            chord_length, setup.airfoil_settings, "Chord", None
        )
        if self._chord is None:
            raise ValueError(
                "Chord length must be provided either directly or in the "
                "setup."
            )

        if altitude is not None and (
                temperature is None or pressure is None or density is None):
            SimpleLogger.warning(
                "Altitude provided, calculating temperature, pressure, and "
                "density using ISA if missing."
            )

            convert_to_kinematic = False
            if pressure is None and self._kinematic_pressure:
                convert_to_kinematic = True

            temperature, pressure, density = (
                ph.international_standard_atmosphere(
                    altitude, temperature, pressure, density)
            )

            SimpleLogger.log(
                f"Conditions: {altitude} m, temperature: {temperature} K, "
                f"absolute pressure: {pressure} Pa, ρ={density} kg/m³"
            )

            if convert_to_kinematic:
                pressure *= density
                SimpleLogger.warning(
                    "Converting calculated absolute pressure to kinematic "
                    "pressure."
                )

            self._temperature = temperature
            self._pressure = pressure
            self._density = density

        if (
            self.nu is None and self.temperature is not None
            and self.density is not None
        ):
            SimpleLogger.warning(
                "Calculating kinematic viscosity based on temperature and "
                "density."
            )

            self._nu = ph.kinematic_viscosity_air(
                self.temperature, self.density)

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
            self.from_mach_number(
                self.mach_number, self.pressure, self.temperature)
        elif (
            self.reynolds_number is not None and self.nu is not None
            and self.pressure is not None
        ):
            self._as_reynolds = True
            self.from_reynolds_number(
                self.reynolds_number, self.nu, self.pressure)
        else:
            SimpleLogger.warning(
                "Insufficient parameters to determine boundary conditions. No "
                "parameters initialized automatically."
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
            "Velocity provided, setting boundary conditions from it "
            "(default variant)."
        )
        if velocity is not None:
            self._velocity = velocity
        if pressure is not None:
            self._pressure = pressure

        reynolds_number = (np.linalg.norm(self.velocity) * self.chord) / self.nu
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
            "Mach number provided, calculating velocity and setting boundary "
            "conditions."
        )
        if mach_number is not None:
            self._mach_number = mach_number
        if pressure is not None:
            self._pressure = pressure
        if temperature is not None:
            self._temperature = temperature

        speed_of_sound = np.sqrt(1.4 * 287.05 * self.temperature)
        self._velocity = np.array([mach_number * speed_of_sound, 0.0, 0.0])

        reynolds_number = (np.linalg.norm(self.velocity) * self.chord) / self.nu
        self._reynolds_number = reynolds_number

    def from_reynolds_number(
            self,
            reynolds_number: float = None,
            nu: float = None,
            pressure: float = None,
    ) -> None:
        """
        Calculate velocity from Reynolds number if it is not provided directly
        and set boundary conditions

        Args:
            reynolds_number (float): Reynolds number.
            nu (float): Kinematic viscosity.
            pressure (float): Pressure value.
        """
        SimpleLogger.log(
            "Reynolds number provided, calculating velocity and setting "
            "boundary conditions."
        )
        if reynolds_number is not None:
            self._reynolds_number = reynolds_number
        if nu is not None:
            self._nu = nu
        if pressure is not None:
            self._pressure = pressure

        self._velocity = np.array([self.reynolds_number * self.nu / self.chord, 0.0, 0.0])

        if self.temperature is not None:
            speed_of_sound = np.sqrt(1.4 * 287.05 * self.temperature)
            mach_number = np.linalg.norm(self.velocity) / speed_of_sound
            self._mach_number = mach_number
        else:
            self._mach_number = None

    def velocity_bc(self) -> str:
        """
        Generate the formatted string for the '0/U' boundary condition file.
        """
        velocity = self.velocity
        if isinstance(velocity, (list, tuple, np.ndarray)):
            velocity_str = " ".join(str(v) for v in velocity)
        else:
            velocity_str = f"{velocity} 0 0"
        return U_bc(
            velocity_str=velocity_str,
            setup=self.bc_setup
        )

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

    def pressure_bc(self) -> str:
        """
        Generate the formatted string for the '0/p' boundary condition file.
        """
        if self._kinematic_pressure:
            pressure_value = self.pressure
        else:
            SimpleLogger.warning(
                "Pressure boundary condition set as absolute pressure, "
                "converting to kinematic for OpenFOAM."
            )
            pressure_value = self.pressure / self.density

        return p_bc(
            pressure_value=str(pressure_value),
            setup=self.bc_setup
        )

    def turbulence_bc(self) -> dict[str, str]:
        """
        Generate the formatted string for the turbulence boundary conditions
        files for the current turbulence model. Supported models: kOmegaSST,
        kOmegaSSTLM, kEpsilon and SpalartAllmaras.

        Returns:
            dict[str, str]: Formatted strings for the turbulence boundary
                condition files.
        """
        model = self._turbulence_model.lower()

        k_inlet_type = self.bc_setup.get("InletK", "fixedValue") or "fixedValue"
        omega_inlet_type = self.bc_setup.get("InletOmega", "fixedValue") or "fixedValue"

        airfoil_nut = self.bc_setup.get(
            "AirfoilNut", "nutUSpaldingWallFunction"
        ) or "nutUSpaldingWallFunction"

        if model == "komegasst":
            return {
                "k": k_bc(
                    k_value=str(self.k),
                    inlet_type=k_inlet_type,
                    intensity=self._turbulence_intensity,
                    setup=self.bc_setup
                ),
                "omega": omega_bc(
                    omega_value=str(self.omega),
                    inlet_type=omega_inlet_type,
                    mixing_length=self.turbulence_length_scale,
                    setup=self.bc_setup
                ),
                "nut": nut_bc(
                    airfoil_wall_function=airfoil_nut,
                    setup=self.bc_setup
                )
            }
        elif model == "komegasstlm":
            return {
                "k": k_bc(
                    k_value=str(self.k),
                    inlet_type=k_inlet_type,
                    intensity=self._turbulence_intensity,
                    setup=self.bc_setup
                ),
                "omega": omega_bc(
                    omega_value=str(self.omega),
                    inlet_type=omega_inlet_type,
                    mixing_length=self.turbulence_length_scale,
                    setup=self.bc_setup
                ),
                "nut": nut_bc(
                    airfoil_wall_function=airfoil_nut,
                    setup=self.bc_setup
                ),
                "gammaInt": gammaInt_bc(
                    gammaInt="1e-5",
                    setup=self.bc_setup
                ),
                "ReThetat": retheta_bc(
                    retheta="1000",
                    setup=self.bc_setup
                )
            }
        elif model == "kepsilon":
            return {
                "k": k_bc(
                    k_value=str(self.k),
                    inlet_type=k_inlet_type,
                    intensity=self._turbulence_intensity,
                    setup=self.bc_setup
                ),
                "epsilon": epsilon_bc(
                    epsilon=str(self.epsilon),
                    setup=self.bc_setup
                ),
                "nut": nut_bc(
                    airfoil_wall_function=airfoil_nut,
                    setup=self.bc_setup
                )
            }
        elif model == "spalartallmaras":
            return {
                "nut": nut_bc(
                    airfoil_wall_function=airfoil_nut,
                    setup=self.bc_setup
                )
            }
        else:
            raise ValueError(f"Unsupported turbulence model: {model}")

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
            SimpleLogger.log("Boundary conditions defined from velocity")
        elif self._as_mach:
            SimpleLogger.log("Boundary conditions defined from Mach number")
        elif self._as_reynolds:
            SimpleLogger.log(
                "Boundary conditions defined from Reynolds number"
            )

        SimpleLogger.log(
            "Main parameters: \n"
            f"                  Velocity: {self.velocity} m/s\n"
        )
        if self._kinematic_pressure:
            SimpleLogger.log(
                f"                  Pressure: {self.pressure} (kinematic) m²/s²"
            )
        else:
            SimpleLogger.log(
                f"                  Pressure: {self.pressure} Pa"
            )

        SimpleLogger.log(
            "Additional parameters: \n"
            f"                  Mach Number: {self.mach_number}\n"
            f"                  Reynolds Number: {self.reynolds_number}\n"
            f"                  Altitude: {self.altitude} m\n"
        )
        SimpleLogger.log(
            "Fluid parameters: \n"
            f"                  Kinematic Viscosity: {self.nu} m²/s\n"
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
        if self._as_velocity:
            definition_method = "Velocity"
        elif self._as_mach:
            definition_method = "Mach Number"
        elif self._as_reynolds:
            definition_method = "Reynolds Number"
        else:
            definition_method = "Unknown"

        if isinstance(self.velocity, np.ndarray):
            velocity_magnitude = np.linalg.norm(self.velocity)
        else:
            velocity_magnitude = self.velocity

        pressure_unit = "m²/s²" if self._kinematic_pressure else "Pa"
        pressure_value = self.pressure if self.pressure is not None else "N/A"

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

        data = [
            ["Parameter", "Unit", "Value"],
            ["Definition Method", "-", definition_method],
            ["Velocity", "m/s", velocity_magnitude],
            ["Pressure", pressure_unit, pressure_value],
            ["Mach Number", "-",
             self.mach_number if self.mach_number is not None else "N/A"],
            ["Reynolds Number", "-",
             self.reynolds_number if self.reynolds_number is not None else "N/A"],
            ["Altitude", "m",
             self.altitude if self.altitude is not None else "N/A"],
            ["Kinematic Viscosity", "m2/s",
             self.nu if self.nu is not None else "N/A"],
            ["Density", "kg/m3",
             self.density if self.density is not None else "N/A"],
            ["Temperature", "K",
             self.temperature if self.temperature is not None else "N/A"],
            ["Chord Length", "m", self.chord],
            ["Turbulence Model", "-", self.turbulence_model],
            ["Turbulence Intensity", "-", self.turbulence_intensity],
            ["Turbulence Length Scale", "m", self.turbulence_length_scale],
            ["k", "m2/s2", self.k],
            [turb_param, turb_unit, turb_value]
        ]

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
        if isinstance(self._velocity, (int, float)):
            return np.array([self._velocity, 0.0, 0.0])
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
        Get the turbulent kinetic energy.

        Returns:
            float: The turbulent kinetic energy.
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
        C_mu = 0.09
        omega = np.sqrt(self.k) / (C_mu ** 0.25
                                   * self.turbulence_length_scale)
        return omega

    @property
    def epsilon(self) -> float:
        """
        Get the dissipation rate for k-epsilon model.

        Returns:
            float: The dissipation rate.
        """
        C_mu = 0.09
        epsilon = (C_mu ** 0.75) * (self.k ** 1.5) / (
            self.turbulence_length_scale)
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
        return self._turbulence_length_scale * self.chord
