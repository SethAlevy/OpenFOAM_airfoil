# Code Documentation

Below you can find the generated API documentation grouped by functional areas of the project.

## Airfoil

Implementation of the airfoil generation and exporting utilities. Able to handle NACA4, NACA5 and UIUC downloading.

::: airfoil.airfoil
::: airfoil.airfoil_base

## Initial Settings

Initial JSON configuration parser and validator for simulation setup files.

::: settings.initial_settings

## Simulation

### Case Preparation

Utilities for creating OpenFOAM case directory structure and preparing simulation files.

::: simulation.preparation.case_structure
::: simulation.preparation.prepare_case

### OpenFOAM Configuration

Modules for generating OpenFOAM configuration files, boundary conditions, and mesh dictionaries.

#### Mesh Generation

::: simulation.openfoam.block_mesh
::: simulation.openfoam.cf_mesh_2d
::: simulation.openfoam.snappy_hex_mesh

#### Boundary Conditions and Dictionaries

::: simulation.openfoam.boundary_condition
::: simulation.openfoam.constant_dir
::: simulation.openfoam.system_dir

## Postprocessing

Tools for extracting, aggregating, and visualizing simulation results.

### Summary Generation

::: postprocess.generate_summary
::: postprocess.aggregate_summaries

### Plotting and Visualization

#### Single Case Plots

::: postprocess.plotting.generate_single_case_plots

#### Polar Plots (Multi-case)

::: postprocess.plotting.generate_polar_plots

#### Matplotlib Utilities

::: postprocess.plotting.matplotlib_plots

#### PyVista Utilities

::: postprocess.plotting.pyvista_plots

## Templates

### OpenFOAM Template Files

Pre-configured templates for OpenFOAM dictionary files.

::: templates.openfoam_template_files.block_mesh_files
::: templates.openfoam_template_files.boundary_files
::: templates.openfoam_template_files.cf_mesh_files
::: templates.openfoam_template_files.constant_files
::: templates.openfoam_template_files.snappy_hex_mesh_file
::: templates.openfoam_template_files.system_files

### Python Template Files

Data structures and configuration templates for Python modules.

::: templates.python_template_files.airfoil_template
::: templates.python_template_files.boundary_condition_template
::: templates.python_template_files.initial_settings_template
::: templates.python_template_files.plot_config

## Utils

General utility functions for geometry, physics calculations, logging, and common operations.

::: utils.geometry
::: utils.physics
::: utils.logger
::: utils.utilities


