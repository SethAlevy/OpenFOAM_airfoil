from pathlib import Path
from templates.initial_settings_template import Settings


def transport_properties_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the transportProperties file for OpenFOAM simulation.

    Args:
        setup (Settings): The simulation settings.
        output_path (Path): The path to save the transportProperties file.
    """
    fluid = setup.simulation_settings.get("Fluid", {})
    transport_model = fluid.get("TransportModel", "Newtonian")
    nu = fluid.get("KinematicViscosity", 1.5e-5)
    content = generate_transport_properties_dict(transport_model=transport_model, nu=nu)
    with open(output_path, "w") as file:
        file.write(content)


def generate_transport_properties_dict(
        transport_model: str = "Newtonian",
        nu: float = 1.5e-5
) -> str:
    """
    Generate transportProperties file content.

    Args:
        transport_model (str): The transport model.
        nu (float): The kinematic viscosity.

    Returns:
        str: The filled transportProperties content.
    """
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}}\n
transportModel  {transport_model};\n
nu              [0 2 -1 0 0 0 0] {nu};"""


def model_specific_coeffs_str(setup: Settings, model: str) -> str:
    """
    Generate a string for model-specific coefficients for the turbulenceProperties file.

    Args:
        setup (Settings): The simulation settings.
        model (str): The turbulence model name (e.g., "kOmegaSST", "kEpsilon").

    Returns:
        str: The formatted string for model-specific coefficients, or empty if none.
    """
    turb = setup.simulation_settings.get("Turbulence", {})
    coeffs = turb.get(model, {})
    if not coeffs:
        return ""
    coeffs_block_name = f"{model}Coeffs"
    coeffs_lines = []
    for k, v in coeffs.items():
        coeffs_lines.append(f"        {k}    {v};")
    coeffs_str = "\n".join(coeffs_lines)
    return f"\n    {coeffs_block_name}\n    {{\n{coeffs_str}\n    }}"


def generate_turbulence_properties_dict(
        simulation_type: str = "RAS",
        turbulence: str = "on",
        print_coeffs: str = "on",
        RAS_model: str = "kOmegaSST",
        kappa: float = 0.41,
        E: float = 9.8,
        model_coeffs: str = ""
) -> str:
    """
    Generate turbulenceProperties file content.

    Args:
        simulation_type (str): The simulation type.
        turbulence (str): Turbulence setting.
        print_coeffs (str): Print coefficients setting.
        RAS_model (str): The RAS turbulence model.
        kappa (float): von Kármán constant.
        E (float): log-law constant.
        model_coeffs (str): Model-specific coefficients block.

    Returns:
        str: The filled turbulenceProperties content.
    """
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}}\n
simulationType  {simulation_type};\n
{simulation_type}\n
{{
    turbulence      {turbulence};
    printCoeffs     {print_coeffs};
    RASModel        {RAS_model};
    kappa           {kappa};
    E               {E};{model_coeffs}
}}"""


def turbulence_properties_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the turbulenceProperties file for OpenFOAM simulation based on the provided
    settings.

    Args:
        setup (Settings): The simulation settings.
        output_path (Path): The path to save the turbulenceProperties file.
    """
    turbulence_settings = setup.simulation_settings.get("Turbulence", {})
    simulation_type = turbulence_settings.get("SimulationType", "RAS")
    turbulence = turbulence_settings.get("Turbulence", "on")
    print_coeffs = turbulence_settings.get("PrintCoeffs", "on")
    model = turbulence_settings.get("Model", "kOmegaSST")
    kappa = turbulence_settings.get("Kappa", 0.41)
    E = turbulence_settings.get("E", 9.8)
    model_coeffs = model_specific_coeffs_str(setup, model)
    content = generate_turbulence_properties_dict(
        simulation_type=simulation_type,
        turbulence=turbulence,
        print_coeffs=print_coeffs,
        RAS_model=model,
        kappa=kappa,
        E=E,
        model_coeffs=model_coeffs
    )
    with open(output_path, "w") as file:
        file.write(content)
