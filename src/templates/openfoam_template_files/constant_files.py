"""
OpenFOAM constant directory dictionary template generation.

This module contains functions to generate OpenFOAM constant directory files:
- transportProperties
- turbulenceProperties
"""


def generate_transport_properties_dict(
    transport_model: str = "Newtonian",
    nu: float = 1.5e-5
) -> str:
    """
    Generate transportProperties file content.

    Defines fluid transport properties for incompressible flow simulations.

    Args:
        transport_model (str): Transport model type (e.g., "Newtonian").
        nu (float): Kinematic viscosity [m²/s].

    Returns:
        str: Complete transportProperties file content.
    """
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}}

transportModel  {transport_model};

nu              [0 2 -1 0 0 0 0] {nu};
"""


def model_coeffs_section(coeffs_dict: dict, model_name: str) -> str:
    """
    Generate model-specific coefficients section.

    Args:
        coeffs_dict (dict): Dictionary of coefficient names and values.
        model_name (str): Name of the turbulence model (e.g., "kOmegaSST").

    Returns:
        str: Formatted coefficients block, or empty string if no coeffs.
    """
    if not coeffs_dict:
        return ""

    coeffs_lines = [f"        {k}    {v};" for k, v in coeffs_dict.items()]
    coeffs_str = "\n".join(coeffs_lines)
    coeffs_block_name = f"{model_name}Coeffs"

    return f"""
    {coeffs_block_name}
    {{
{coeffs_str}
    }}"""


def generate_turbulence_properties_dict(
    simulation_type: str = "RAS",
    turbulence: str = "on",
    print_coeffs: str = "on",
    ras_model: str = "kOmegaSST",
    kappa: float = 0.41,
    E: float = 9.8,
    model_coeffs: str = ""
) -> str:
    """
    Generate turbulenceProperties file content.

    Defines turbulence modeling settings for RANS/LES simulations.

    Args:
        simulation_type (str): Simulation type ("RAS", "LES", or "DNS").
        turbulence (str): Enable turbulence modeling ("on" or "off").
        print_coeffs (str): Print model coefficients ("on" or "off").
        ras_model (str): RANS turbulence model name.
        kappa (float): von Kármán constant (typically 0.41).
        E (float): Wall law constant (typically 9.8).
        model_coeffs (str): Model-specific coefficients block (pre-formatted).

    Returns:
        str: Complete turbulenceProperties file content.
    """
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}}

simulationType  {simulation_type};

{simulation_type}
{{
    turbulence      {turbulence};
    printCoeffs     {print_coeffs};
    RASModel        {ras_model};
    kappa           {kappa};
    E               {E};{model_coeffs}
}}
"""
