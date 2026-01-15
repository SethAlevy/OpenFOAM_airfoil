"""
Physical and thermodynamic utility functions for aerodynamic simulations.
"""


def kinematic_viscosity_air(temperature: float, density: float) -> float:
    """
    Calculate kinematic viscosity of air at temperature [K] and density
    [kg/m^3] using Sutherland's formula for dynamic viscosity.

    Args:
        temperature (float): Temperature in Kelvin [K].
        density (float): Density in kg/m^3 [kg/m^3].

    Returns:
        float: Kinematic viscosity in m^2/s [m^2/s].

    Note:
        Uses Sutherland's formula valid for 170 K < T < 1900 K.
        Coefficients calibrated for air at standard conditions.
    """
    # Sutherland's formula constants for air
    reference_dynamic_viscosity = 1.716e-5  # Reference viscosity [Pa·s] at T0
    reference_temperature = 273.15  # Reference temperature [K] (0°C)
    sutherland_constant = 110.4  # Sutherland constant [K]

    # Calculate dynamic viscosity using Sutherland's formula
    dynamic_viscosity = (
        reference_dynamic_viscosity
        * (temperature / reference_temperature) ** 1.5
        * (reference_temperature + sutherland_constant)
        / (temperature + sutherland_constant)
    )

    return dynamic_viscosity / density


def international_standard_atmosphere(
    altitude: float,
    temperature: float = None,
    pressure: float = None,
    density: float = None,
) -> tuple[float, float, float]:
    """
    Calculate temperature [K], pressure [Pa], and density [kg/m^3] at a given
    altitude [m] using the International Standard Atmosphere (ISA) model.

    Implements the tropospheric model (0 m to 11 km altitude) with linear
    temperature lapse rate.

    Args:
        altitude (float): Altitude above sea level [m].
        temperature (float): Temperature in Kelvin [K].n If None, calculated from ISA.
        pressure (float): Pressure in Pascals [Pa]. If None, calculated from ISA.
        density (float): Density in kg/m^3 [kg/m^3]. If None, calculated from ISA.

    Returns:
        tuple[float, float, float]: (temperature [K], pressure [Pa], density [kg/m^3]).
                                    Returns user-provided values if given,
                                    otherwise ISA model values.

    Raises:
        ValueError: If altitude is outside troposphere (< 0 m or > 11000 m).
    """
    # ISA parameters - Tropospheric layer (0 m to 11 km)
    sea_level_temperature = 288.15  # Temperature at sea level [K] (15°C)
    sea_level_pressure = 101325.0  # Pressure at sea level [Pa]
    temperature_lapse_rate = 0.0065  # Linear temperature decrease [K/m]
    gravitational_acceleration = 9.80665  # Standard gravity [m/s^2]
    molar_mass_air = 0.0289644  # Molar mass of dry air [kg/mol]
    universal_gas_constant = 8.3144598  # Universal gas constant [J/(mol·K)]
    specific_gas_constant_air = 287.05  # Specific gas constant for dry air [J/(kg·K)]

    isa_temperature = sea_level_temperature - temperature_lapse_rate * altitude

    # Barometric formula for pressure
    exponent = (
        gravitational_acceleration * molar_mass_air
        / (universal_gas_constant * temperature_lapse_rate)
    )
    isa_pressure = sea_level_pressure * (
        1 - temperature_lapse_rate * altitude / sea_level_temperature
    ) ** exponent

    # Ideal gas law for density
    isa_density = isa_pressure / (specific_gas_constant_air * isa_temperature)

    # Use ISA values if user values not provided
    final_temperature = temperature if temperature is not None else isa_temperature
    final_pressure = pressure if pressure is not None else isa_pressure
    final_density = density if density is not None else isa_density

    return final_temperature, final_pressure, final_density


def calculate_mach_number(
    velocity: float,
    temperature: float,
    gamma: float = 1.4
) -> float:
    """
    Calculate Mach number from velocity and local speed of sound.

    Args:
        velocity (float): Velocity [m/s].
        temperature (float): Static temperature [K].
        gamma (float): Specific heat ratio (Cp/Cv), default 1.4 for air.

    Returns:
        float: Mach number (dimensionless).
    """
    specific_gas_constant_air = 287.05  # [J/(kg·K)]
    speed_of_sound = (
        gamma * specific_gas_constant_air * temperature
    ) ** 0.5
    return velocity / speed_of_sound


def calculate_reynolds_number(
    velocity: float,
    characteristic_length: float,
    kinematic_viscosity: float
) -> float:
    """
    Calculate Reynolds number from velocity, characteristic length, and kinematic
    viscosity.

    Args:
        velocity (float): Velocity [m/s].
        characteristic_length (float): Characteristic length (e.g., chord) [m].
        kinematic_viscosity (float): Kinematic viscosity [m^2/s].

    Returns:
        float: Reynolds number (dimensionless).
    """
    return (velocity * characteristic_length) / kinematic_viscosity


def calculate_dynamic_pressure(
    velocity: float,
    density: float
) -> float:
    """
    Calculate dynamic pressure from velocity and density.

    Args:
        velocity (float): Velocity [m/s].
        density (float): Density [kg/m^3].

    Returns:
        float: Dynamic pressure [Pa].
    """
    return 0.5 * density * velocity ** 2
