"""
OpenFOAM boundary condition file templates.

Provides functions to generate boundary condition files for various fields.
"""

from typing import Literal


def patch(name: str, setup: dict, default: str) -> str:
    """
    Resolve patch name from parameter, setup dict, or default.

    Args:
        name: Patch name provided as parameter.
        setup: Dictionary with optional patch name defaults.
        default: Default patch name if none provided.

    Returns:
        str: Resolved patch name."""
    return name if name is not None else setup.get(default, default)


def generate_foam_field(
    object_name: str,
    dimensions: str,
    internal_value: str,
    boundary_conditions: dict,
    field_class: str = "volScalarField",
) -> str:
    """
    Generic OpenFOAM field file template generator.

    Args:
        object_name: Object name in FoamFile (e.g., "U", "p", "k")
        dimensions: Dimensional set [kg m s K mol A cd]
        internal_value: Internal field value
        boundary_conditions: Dict mapping patch names to BC specs
        field_class: Field class (volScalarField, volVectorField, etc.)

    Returns:
        Formatted OpenFOAM field file as string
    """
    foam_file_header = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {field_class};
    object      {object_name};
}}

dimensions      {dimensions};

internalField   uniform {internal_value};

boundaryField
{{"""

    boundary_entries = "\n".join(
        f"    {bc_spec['patch_name']}\n    {{\n{bc_spec['content']}\n    }}"
        for bc_spec in boundary_conditions.values()
    )

    foam_file_footer = """
}
"""
    return foam_file_header + "\n" + boundary_entries + foam_file_footer


def build_bc_spec(patch_name: str, bc_type: str, value: str = None) -> dict:
    """
    Build a boundary condition specification.

    Args:
        patch_name: Name of the patch
        bc_type: OpenFOAM BC type (fixedValue, zeroGradient, etc.)
        value: Value line if needed (e.g., "value   uniform 0.5;")

    Returns:
        Dictionary with patch_name and content keys
    """
    content = f"        type            {bc_type};"
    if value:
        content += f"\n        {value}"

    return {
        "patch_name": patch_name,
        "content": content,
    }


def U_bc(
    velocity_str: str,
    inlet_patch: str = None,
    outlet_patch: str = None,
    lower_wall_patch: str = None,
    upper_wall_patch: str = None,
    front_patch: str = None,
    back_patch: str = None,
    airfoil_patch: str = None,
    setup: dict = None,
) -> str:
    """
    Generate velocity (U) boundary condition file.

    Args:
        velocity_str: Velocity vector string (e.g., "10 0 0")
        inlet_patch: Inlet patch name override
        outlet_patch: Outlet patch name override
        lower_wall_patch: Lower wall patch name override
        upper_wall_patch: Upper wall patch name override
        front_patch: Front patch name override
        back_patch: Back patch name override
        airfoil_patch: Airfoil patch name override
        setup: dict = None,

    Returns:
        OpenFOAM U field file content
    """
    setup = setup or {}
    inlet = patch(inlet_patch, setup, "inlet")
    outlet = patch(outlet_patch, setup, "outlet")
    lower_wall = patch(lower_wall_patch, setup, "lowerWall")
    upper_wall = patch(upper_wall_patch, setup, "upperWall")
    front = patch(front_patch, setup, "front")
    back = patch(back_patch, setup, "back")
    airfoil = patch(airfoil_patch, setup, "airfoil")

    boundary_conditions = {
        "inlet": build_bc_spec(
            inlet, "fixedValue", f"value           uniform ({velocity_str});"
        ),
        "outlet": build_bc_spec(outlet, "zeroGradient"),
        "lowerWall": build_bc_spec(lower_wall, "zeroGradient"),
        "upperWall": build_bc_spec(upper_wall, "zeroGradient"),
        "front": build_bc_spec(front, "symmetryPlane"),
        "back": build_bc_spec(back, "symmetryPlane"),
        "airfoil": build_bc_spec(airfoil, "noSlip"),
    }

    return generate_foam_field(
        "U",
        "[0 1 -1 0 0 0 0]",
        f"({velocity_str})",
        boundary_conditions,
        "volVectorField",
    )


def p_bc(
    pressure_value: str,
    inlet_patch : str = None,
    outlet_patch: str = None,
    lower_wall_patch: str = None,
    upper_wall_patch: str = None,
    front_patch: str = None,
    back_patch: str = None,
    airfoil_patch: str = None,
    setup: dict = None,
) -> str:
    """
    Generate pressure (p) boundary condition file.

    Args:
        pressure_value: Pressure value string
        inlet_patch: Inlet patch name override
        outlet_patch: Outlet patch name override
        lower_wall_patch: Lower wall patch name override
        upper_wall_patch: Upper wall patch name override
        front_patch: Front patch name override
        back_patch: Back patch name override
        airfoil_patch: Airfoil patch name override
        setup: Dictionary with optional patch name defaults

    Returns:
        OpenFOAM p field file content
    """
    setup = setup or {}
    inlet = patch(inlet_patch, setup, "inlet")
    outlet = patch(outlet_patch, setup, "outlet")
    lower_wall = patch(lower_wall_patch, setup, "lowerWall")
    upper_wall = patch(upper_wall_patch, setup, "upperWall")
    front = patch(front_patch, setup, "front")
    back = patch(back_patch, setup, "back")
    airfoil = patch(airfoil_patch, setup, "airfoil")

    boundary_conditions = {
        "inlet": build_bc_spec(inlet, "zeroGradient"),
        "outlet": build_bc_spec(
            outlet, "fixedValue", f"value           uniform {pressure_value};"
        ),
        "lowerWall": build_bc_spec(lower_wall, "zeroGradient"),
        "upperWall": build_bc_spec(upper_wall, "zeroGradient"),
        "front": build_bc_spec(front, "symmetryPlane"),
        "back": build_bc_spec(back, "symmetryPlane"),
        "airfoil": build_bc_spec(airfoil, "zeroGradient"),
    }

    return generate_foam_field(
        "p", "[0 2 -2 0 0 0 0]", pressure_value, boundary_conditions
    )


def k_bc(
    k_value: str,
    inlet_type: Literal[
        "fixedValue", "turbulentIntensityKineticEnergyInlet"
    ] = "turbulentIntensityKineticEnergyInlet",
    intensity: float = None,
    inlet_patch: str = None,
    outlet_patch: str = None,
    lower_wall_patch: str = None,
    upper_wall_patch: str = None,
    front_patch: str = None,
    back_patch: str = None,
    airfoil_patch: str = None,
    setup: dict = None,
) -> str:
    """
    Generate turbulent kinetic energy (k) boundary condition file.

    Args:
        k_value: Pre-calculated k value string
        inlet_type: Inlet BC type ('fixedValue' or
                   'turbulentIntensityKineticEnergyInlet')
        intensity: Turbulence intensity (required if inlet_type is
                  'turbulentIntensityKineticEnergyInlet')
        inlet_patch: Inlet patch name override
        outlet_patch: Outlet patch name override
        lower_wall_patch: Lower wall patch name override
        upper_wall_patch: Upper wall patch name override
        front_patch: Front patch name override
        back_patch: Back patch name override
        airfoil_patch: Airfoil patch name override
        setup: Dictionary with optional patch name defaults

    Returns:
        OpenFOAM k field file content

    Raises:
        ValueError: If inlet_type requires intensity but not provided
    """
    setup = setup or {}
    inlet = patch(inlet_patch, setup, "inlet")
    outlet = patch(outlet_patch, setup, "outlet")
    lower_wall = patch(lower_wall_patch, setup, "lowerWall")
    upper_wall = patch(upper_wall_patch, setup, "upperWall")
    front = patch(front_patch, setup, "front")
    back = patch(back_patch, setup, "back")
    airfoil = patch(airfoil_patch, setup, "airfoil")

    # Handle empty string inlet_type
    if not inlet_type or inlet_type == "fixedValue":
        inlet_bc = build_bc_spec(
            inlet, "fixedValue", f"value           uniform {k_value};"
        )
    elif inlet_type == "turbulentIntensityKineticEnergyInlet":
        if intensity is None:
            raise ValueError(
                "intensity required for "
                "inlet_type='turbulentIntensityKineticEnergyInlet'"
            )
        inlet_bc = build_bc_spec(
            inlet,
            "turbulentIntensityKineticEnergyInlet",
            f"intensity       {intensity};\n"
            f"        value           uniform {k_value};",
        )
    else:
        raise ValueError(f"Unsupported inlet_type for k: {inlet_type}")

    boundary_conditions = {
        "inlet": inlet_bc,
        "outlet": build_bc_spec(outlet, "zeroGradient"),
        "lowerWall": build_bc_spec(lower_wall, "zeroGradient"),
        "upperWall": build_bc_spec(upper_wall, "zeroGradient"),
        "front": build_bc_spec(front, "symmetryPlane"),
        "back": build_bc_spec(back, "symmetryPlane"),
        "airfoil": build_bc_spec(
            airfoil, "kqRWallFunction", "value           uniform 0;"
        ),
    }

    return generate_foam_field(
        "k", "[0 2 -2 0 0 0 0]", k_value, boundary_conditions
    )


def omega_bc(
    omega_value: str,
    inlet_type: Literal[
        "fixedValue", "turbulentMixingLengthFrequencyInlet"
    ] = "turbulentMixingLengthFrequencyInlet",
    mixing_length: float = None,
    inlet_patch: str = None,
    outlet_patch: str = None,
    lower_wall_patch: str = None,
    upper_wall_patch: str = None,
    front_patch: str = None,
    back_patch: str = None,
    airfoil_patch: str = None,
    setup: dict = None,
) -> str:
    """
    Generate specific dissipation rate (omega) boundary condition file.

    Args:
        omega_value: Pre-calculated omega value string
        inlet_type: Inlet BC type ('fixedValue' or
                   'turbulentMixingLengthFrequencyInlet')
        mixing_length: Turbulent mixing length (required if inlet_type is
                      'turbulentMixingLengthFrequencyInlet')
        inlet_patch: Inlet patch name override
        outlet_patch: Outlet patch name override
        lower_wall_patch: Lower wall patch name override
        upper_wall_patch: Upper wall patch name override
        front_patch: Front patch name override
        back_patch: Back patch name override
        airfoil_patch: Airfoil patch name override
        setup: Dictionary with optional patch name defaults

    Returns:
        OpenFOAM omega field file content

    Raises:
        ValueError: If inlet_type requires mixing_length but not provided
    """
    setup = setup or {}
    inlet = patch(inlet_patch, setup, "inlet")
    outlet = patch(outlet_patch, setup, "outlet")
    lower_wall = patch(lower_wall_patch, setup, "lowerWall")
    upper_wall = patch(upper_wall_patch, setup, "upperWall")
    front = patch(front_patch, setup, "front")
    back = patch(back_patch, setup, "back")
    airfoil = patch(airfoil_patch, setup, "airfoil")

    # Handle empty string inlet_type
    if not inlet_type or inlet_type == "fixedValue":
        inlet_bc = build_bc_spec(
            inlet, "fixedValue", f"value           uniform {omega_value};"
        )
    elif inlet_type == "turbulentMixingLengthFrequencyInlet":
        if mixing_length is None:
            raise ValueError(
                "mixing_length required for "
                "inlet_type='turbulentMixingLengthFrequencyInlet'"
            )
        inlet_bc = build_bc_spec(
            inlet,
            "turbulentMixingLengthFrequencyInlet",
            f"mixingLength    {mixing_length};\n"
            f"        value           uniform {omega_value};",
        )
    else:
        raise ValueError(f"Unsupported inlet_type for omega: {inlet_type}")

    boundary_conditions = {
        "inlet": inlet_bc,
        "outlet": build_bc_spec(outlet, "zeroGradient"),
        "lowerWall": build_bc_spec(lower_wall, "zeroGradient"),
        "upperWall": build_bc_spec(upper_wall, "zeroGradient"),
        "front": build_bc_spec(front, "symmetryPlane"),
        "back": build_bc_spec(back, "symmetryPlane"),
        "airfoil": build_bc_spec(
            airfoil, "omegaWallFunction", "value           uniform 0;"
        ),
    }

    return generate_foam_field(
        "omega", "[0 0 -1 0 0 0 0]", omega_value, boundary_conditions
    )


def epsilon_bc(
    epsilon: str,
    inlet_patch: str = None,
    outlet_patch: str = None,
    lower_wall_patch: str = None,
    upper_wall_patch: str = None,
    front_patch: str = None,
    back_patch: str = None,
    airfoil_patch: str = None,
    setup: dict = None,
) -> str:
    """
    Generate turbulent dissipation rate (epsilon) boundary condition file.

    Args:
        epsilon: Epsilon value string
        inlet_patch: Inlet patch name override
        outlet_patch: Outlet patch name override
        lower_wall_patch: Lower wall patch name override
        upper_wall_patch: Upper wall patch name override
        front_patch: Front patch name override
        back_patch: Back patch name override
        airfoil_patch: Airfoil patch name override
        setup: Dictionary with optional patch name defaults

    Returns:
        OpenFOAM epsilon field file content
    """
    setup = setup or {}
    inlet = patch(inlet_patch, setup, "inlet")
    outlet = patch(outlet_patch, setup, "outlet")
    lower_wall = patch(lower_wall_patch, setup, "lowerWall")
    upper_wall = patch(upper_wall_patch, setup, "upperWall")
    front = patch(front_patch, setup, "front")
    back = patch(back_patch, setup, "back")
    airfoil = patch(airfoil_patch, setup, "airfoil")

    boundary_conditions = {
        "inlet": build_bc_spec(
            inlet, "fixedValue", f"value           uniform {epsilon};"
        ),
        "outlet": build_bc_spec(outlet, "zeroGradient"),
        "lowerWall": build_bc_spec(lower_wall, "zeroGradient"),
        "upperWall": build_bc_spec(upper_wall, "zeroGradient"),
        "front": build_bc_spec(
            front, "symmetryPlane", f"value           uniform {epsilon};"
        ),
        "back": build_bc_spec(
            back, "symmetryPlane", f"value           uniform {epsilon};"
        ),
        "airfoil": build_bc_spec(
            airfoil, "fixedValue", "value           uniform 0;"
        ),
    }

    return generate_foam_field(
        "epsilon", "[0 2 -3 0 0 0 0]", epsilon, boundary_conditions
    )


def nut_bc(
    airfoil_wall_function: Literal[
        "nutUSpaldingWallFunction", "nutLowReWallFunction"
    ] = "nutLowReWallFunction",
    inlet_patch: str = None,
    outlet_patch: str = None,
    lower_wall_patch: str = None,
    upper_wall_patch: str = None,
    front_patch: str = None,
    back_patch: str = None,
    airfoil_patch: str = None,
    setup: dict = None,
) -> str:
    """
    Generate turbulent viscosity (nut) boundary condition file.

    Args:
        airfoil_wall_function: Wall function for airfoil ('nutUSpaldingWallFunction'
                              or 'nutLowReWallFunction').
        inlet_patch: Inlet patch name override
        outlet_patch: Outlet patch name override
        lower_wall_patch: Lower wall patch name override
        upper_wall_patch: Upper wall patch name override
        front_patch: Front patch name override
        back_patch: Back patch name override
        airfoil_patch: Airfoil patch name override
        setup: Dictionary with optional patch name defaults

    Returns:
        OpenFOAM nut field file content

    Raises:
        ValueError: If airfoil_wall_function is not a supported wall function
    """
    setup = setup or {}
    inlet = patch(inlet_patch, setup, "inlet")
    outlet = patch(outlet_patch, setup, "outlet")
    lower_wall = patch(lower_wall_patch, setup, "lowerWall")
    upper_wall = patch(upper_wall_patch, setup, "upperWall")
    front = patch(front_patch, setup, "front")
    back = patch(back_patch, setup, "back")
    airfoil = patch(airfoil_patch, setup, "airfoil")

    supported_wall_functions = [
        "nutUSpaldingWallFunction",
        "nutLowReWallFunction"
    ]
    if airfoil_wall_function not in supported_wall_functions:
        raise ValueError(
            f"Unsupported wall function: {airfoil_wall_function}. "
            f"Supported: {', '.join(supported_wall_functions)}"
        )

    boundary_conditions = {
        "inlet": build_bc_spec(inlet, "zeroGradient"),
        "outlet": build_bc_spec(outlet, "zeroGradient"),
        "lowerWall": build_bc_spec(lower_wall, "zeroGradient"),
        "upperWall": build_bc_spec(upper_wall, "zeroGradient"),
        "front": build_bc_spec(front, "symmetryPlane"),
        "back": build_bc_spec(back, "symmetryPlane"),
        "airfoil": build_bc_spec(
            airfoil, airfoil_wall_function, "value           uniform 0;"
        ),
    }

    return generate_foam_field(
        "nut", "[0 2 -1 0 0 0 0]", "0", boundary_conditions
    )


def gammaInt_bc(
    gammaInt: str,
    inlet_patch: str = None,
    outlet_patch: str = None,
    lower_wall_patch: str = None,
    upper_wall_patch: str = None,
    front_patch: str = None,
    back_patch: str = None,
    airfoil_patch: str = None,
    setup: dict = None,
) -> str:
    """
    Generate intermittency (gammaInt) boundary condition file.

    Args:
        gammaInt: GammaInt value string
        inlet_patch: Inlet patch name override
        outlet_patch: Outlet patch name override
        lower_wall_patch: Lower wall patch name override
        upper_wall_patch: Upper wall patch name override
        front_patch: Front patch name override
        back_patch: Back patch name override
        airfoil_patch: Airfoil patch name override
        setup: Dictionary with optional patch name defaults

    Returns:
        OpenFOAM gammaInt field file content
    """
    setup = setup or {}
    inlet = patch(inlet_patch, setup, "inlet")
    outlet = patch(outlet_patch, setup, "outlet")
    lower_wall = patch(lower_wall_patch, setup, "lowerWall")
    upper_wall = patch(upper_wall_patch, setup, "upperWall")
    front = patch(front_patch, setup, "front")
    back = patch(back_patch, setup, "back")
    airfoil = patch(airfoil_patch, setup, "airfoil")

    boundary_conditions = {
        "inlet": build_bc_spec(
            inlet, "fixedValue", f"value           uniform {gammaInt};"
        ),
        "outlet": build_bc_spec(outlet, "zeroGradient"),
        "lowerWall": build_bc_spec(lower_wall, "zeroGradient"),
        "upperWall": build_bc_spec(upper_wall, "zeroGradient"),
        "front": build_bc_spec(front, "symmetryPlane"),
        "back": build_bc_spec(back, "symmetryPlane"),
        "airfoil": build_bc_spec(airfoil, "zeroGradient"),
    }

    return generate_foam_field(
        "gammaInt",
        "[0 0 0 0 0 0 0]",
        gammaInt,
        boundary_conditions,
    )


def retheta_bc(
    retheta: str,
    inlet_patch: str = None,
    outlet_patch: str = None,
    lower_wall_patch: str = None,
    upper_wall_patch: str = None,
    front_patch: str = None,
    back_patch: str = None,
    airfoil_patch: str = None,
    setup: dict = None,
) -> str:
    """
    Generate Reynolds theta (ReTheta) boundary condition file.

    Args:
        retheta: ReTheta value string
        inlet_patch: Inlet patch name override
        outlet_patch: Outlet patch name override
        lower_wall_patch: Lower wall patch name override
        upper_wall_patch: Upper wall patch name override
        front_patch: Front patch name override
        back_patch: Back patch name override
        airfoil_patch: Airfoil patch name override
        setup: Dictionary with optional patch name defaults

    Returns:
        OpenFOAM ReTheta field file content
    """
    setup = setup or {}
    inlet = patch(inlet_patch, setup, "inlet")
    outlet = patch(outlet_patch, setup, "outlet")
    lower_wall = patch(lower_wall_patch, setup, "lowerWall")
    upper_wall = patch(upper_wall_patch, setup, "upperWall")
    front = patch(front_patch, setup, "front")
    back = patch(back_patch, setup, "back")
    airfoil = patch(airfoil_patch, setup, "airfoil")

    boundary_conditions = {
        "inlet": build_bc_spec(
            inlet, "fixedValue", f"value           uniform {retheta};"
        ),
        "outlet": build_bc_spec(outlet, "zeroGradient"),
        "lowerWall": build_bc_spec(lower_wall, "zeroGradient"),
        "upperWall": build_bc_spec(upper_wall, "zeroGradient"),
        "front": build_bc_spec(front, "symmetryPlane"),
        "back": build_bc_spec(back, "symmetryPlane"),
        "airfoil": build_bc_spec(
            airfoil, "fixedValue", "value           uniform 1000;"
        ),
    }

    return generate_foam_field(
        "ReTheta", "[0 0 0 0 0 0 0]", retheta, boundary_conditions
    )
