from pathlib import Path
from templates.initial_settings_template import Settings


def transport_properties_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the transportProperties file for OpenFOAM simulation.
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


def turbulence_properties_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the turbulenceProperties file for OpenFOAM simulation.
    """
    turb = setup.simulation_settings.get("Turbulence", {})
    simulation_type = turb.get("SimulationType", "RAS")
    turbulence = turb.get("Turbulence", "on")
    print_coeffs = turb.get("PrintCoeffs", "on")
    model = turb.get("Model", "kOmegaSST")
    content = generate_turbulence_properties_dict(
        simulation_type=simulation_type,
        turbulence=turbulence,
        print_coeffs=print_coeffs,
        RAS_model=model
    )
    with open(output_path, "w") as file:
        file.write(content)


def generate_turbulence_properties_dict(
        simulation_type: str = "RAS",
        turbulence: str = "on",
        print_coeffs: str = "on",
        RAS_model: str = "kOmegaSST"
) -> str:
    """
    Generate turbulenceProperties file content.
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
}}"""
