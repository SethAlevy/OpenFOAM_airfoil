from pathlib import Path
from templates.python_template_files.initial_settings_template import Settings
from templates.openfoam_template_files.constant_files import (
    generate_transport_properties_dict,
    generate_turbulence_properties_dict,
    model_coeffs_section
)


def transport_properties_dict(setup: Settings, output_path: Path) -> None:
    """
    Generate and write transportProperties file for OpenFOAM simulation.

    Args:
        setup (Settings): The simulation settings.
        output_path (Path): The path to save the transportProperties file.
    """
    fluid = setup.simulation_settings.get("Fluid", {})
    transport_model = fluid.get("TransportModel", "Newtonian")
    nu = fluid.get("KinematicViscosity", 1.5e-5)

    content = generate_transport_properties_dict(
        transport_model=transport_model,
        nu=nu
    )

    with open(output_path, "w") as file:
        file.write(content)


def turbulence_properties_dict(setup: Settings, output_path: Path) -> None:
    """
    Generate and write turbulenceProperties file for OpenFOAM simulation.

    Args:
        setup (Settings): The simulation settings.
        output_path (Path): The path to save the turbulenceProperties file.
    """
    turbulence_settings = setup.simulation_settings.get("TurbulenceProperties", {})
    simulation_type = turbulence_settings.get("SimulationType", "RAS")
    turbulence = turbulence_settings.get("Turbulence", "on")
    print_coeffs = turbulence_settings.get("PrintCoeffs", "on")
    model = turbulence_settings.get("Model", "kOmegaSST")
    kappa = turbulence_settings.get("Kappa", 0.41)
    E = turbulence_settings.get("E", 9.8)

    model_coeffs_dict = turbulence_settings.get(model, {})
    model_coeffs_str = model_coeffs_section(model_coeffs_dict, model)

    content = generate_turbulence_properties_dict(
        simulation_type=simulation_type,
        turbulence=turbulence,
        print_coeffs=print_coeffs,
        ras_model=model,
        kappa=kappa,
        E=E,
        model_coeffs=model_coeffs_str
    )

    with open(output_path, "w") as file:
        file.write(content)
