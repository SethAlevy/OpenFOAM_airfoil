import numpy as np
from pathlib import Path
from templates.initial_settings_template import Settings
from simulation.openfoam.boundary_condition import BoundaryConditions
from templates.openfoam_template_files.system_files import (
    generate_control_dict,
    generate_fv_solution_dict,
    generate_fv_schemes_dict,
    generate_decompose_par_dict,
    generate_surface_feature_extract_dict,
    generate_extrude_mesh_dict,
)


def control_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the controlDict file for OpenFOAM simulation based on the provided
    settings. Then write the file.

    Args:
        setup (Settings): The simulation settings json file.
        output_path (Path): The path to save the controlDict file.
    """
    control_dict_setup = setup.simulation_settings.get("ControlDict", {})
    solver = control_dict_setup.get("Solver", "simpleFoam")
    start_from = control_dict_setup.get("StartFrom", "startTime")
    start_time = control_dict_setup.get("StartTime", 0.0)
    stop_at = control_dict_setup.get("StopAt", "endTime")
    end_time = control_dict_setup.get("EndTime", 1000.0)
    delta_t = control_dict_setup.get("DeltaT", 1.0)
    purge_write = control_dict_setup.get("PurgeWrite", 0)
    write_interval = control_dict_setup.get("WriteInterval", 100)
    write_control = control_dict_setup.get("WriteControl", "timeStep")
    write_format = control_dict_setup.get("WriteFormat", "ascii")
    write_precision = control_dict_setup.get("WritePrecision", 6)
    write_compression = control_dict_setup.get("WriteCompression", "off")
    time_format = control_dict_setup.get("TimeFormat", "general")
    time_precision = control_dict_setup.get("TimePrecision", 6)
    run_time_modifiable = control_dict_setup.get("RunTimeModifiable", "true")

    content = generate_control_dict(
        solver=solver,
        end_time=end_time,
        write_interval=write_interval,
        start_from=start_from,
        start_time=start_time,
        stop_at=stop_at,
        delta_t=delta_t,
        purge_write=purge_write,
        write_control=write_control,
        write_format=write_format,
        write_precision=write_precision,
        write_compression=write_compression,
        time_format=time_format,
        time_precision=time_precision,
        run_time_modifiable=run_time_modifiable
    )
    with open(output_path, "w") as file:
        file.write(content)


def fv_solution_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the fvSolution file for OpenFOAM simulation based on the provided
    settings.

    Args:
        setup (Settings): The simulation settings json file.
        output_path (Path): The path to save the fvSolution file.
    """
    fv_solution_setup = setup.simulation_settings.get("FvSolution", {})

    p_setup = fv_solution_setup.get("p", {})
    p_solver = p_setup.get("Solver", "PCG")
    p_preconditioner = p_setup.get("Preconditioner", "DIC")
    p_smoother = p_setup.get("Smoother", None)
    p_tolerance = p_setup.get("Tolerance", 1e-6)
    p_rel_tol = p_setup.get("RelTol", 0)

    u_setup = fv_solution_setup.get("U", {})
    u_solver = u_setup.get("Solver", "smoothSolver")
    u_preconditioner = u_setup.get("Preconditioner", None)
    u_smoother = u_setup.get("Smoother", "symGaussSeidel")
    u_tolerance = u_setup.get("Tolerance", 1e-5)
    u_rel_tol = u_setup.get("RelTol", 0)

    turb_setup = fv_solution_setup.get("Turbulence", {})
    turb_solver = turb_setup.get("Solver", "smoothSolver")
    turb_preconditioner = turb_setup.get("Preconditioner", None)
    turb_smoother = turb_setup.get("Smoother", "symGaussSeidel")
    turb_tolerance = turb_setup.get("Tolerance", 1e-5)
    turb_rel_tol = turb_setup.get("RelTol", 0)

    re_thetat_setup = fv_solution_setup.get("ReThetat", {})
    rethetat_solver = re_thetat_setup.get("Solver", "smoothSolver")
    rethetat_preconditioner = re_thetat_setup.get("Preconditioner", None)
    rethetat_smoother = re_thetat_setup.get("Smoother", "symGaussSeidel")
    rethetat_tolerance = re_thetat_setup.get("Tolerance", 1e-8)
    rethetat_rel_tol = re_thetat_setup.get("RelTol", 0.1)

    gammaint_setup = fv_solution_setup.get("gammaInt", {})
    gammaint_solver = gammaint_setup.get("Solver", "smoothSolver")
    gammaint_preconditioner = gammaint_setup.get("Preconditioner", None)
    gammaint_smoother = gammaint_setup.get("Smoother", "symGaussSeidel")
    gammaint_tolerance = gammaint_setup.get("Tolerance", 1e-8)
    gammaint_rel_tol = gammaint_setup.get("RelTol", 0.1)

    n_non_ortho_correctors = fv_solution_setup.get(
        "n_non_ortho_correctors", 2)

    relaxation_setup = fv_solution_setup.get("RelaxationFactors", {})
    relaxation_p = relaxation_setup.get("p", 0.3)
    relaxation_U = relaxation_setup.get("U", 0.7)
    relaxation_k = relaxation_setup.get("k", 0.7)
    relaxation_omega = relaxation_setup.get("omega", 0.7)
    relaxation_epsilon = relaxation_setup.get("epsilon", 0.7)
    relaxation_nuTilda = relaxation_setup.get("nuTilda", 0.7)
    relaxation_reThetat = relaxation_setup.get("reThetat", 0.4)
    relaxation_gammaInt = relaxation_setup.get("gammaInt", 0.4)

    cache = fv_solution_setup.get("cache", "grad(U);")

    content = generate_fv_solution_dict(
        p_solver=p_solver,
        p_preconditioner=p_preconditioner,
        p_smoother=p_smoother,
        p_tolerance=p_tolerance,
        p_rel_tol=p_rel_tol,
        U_solver=u_solver,
        U_preconditioner=u_preconditioner,
        U_smoother=u_smoother,
        U_tolerance=u_tolerance,
        U_rel_tol=u_rel_tol,
        turb_solver=turb_solver,
        turb_preconditioner=turb_preconditioner,
        turb_smoother=turb_smoother,
        turb_tolerance=turb_tolerance,
        turb_rel_tol=turb_rel_tol,
        n_non_ortho_correctors=n_non_ortho_correctors,
        relaxation_p=relaxation_p,
        relaxation_U=relaxation_U,
        relaxation_k=relaxation_k,
        relaxation_omega=relaxation_omega,
        relaxation_epsilon=relaxation_epsilon,
        relaxation_nuTilda=relaxation_nuTilda,
        cache=cache,
        rethetat_solver=rethetat_solver,
        rethetat_preconditioner=rethetat_preconditioner,
        rethetat_smoother=rethetat_smoother,
        rethetat_tolerance=rethetat_tolerance,
        rethetat_rel_tol=rethetat_rel_tol,
        gammaint_solver=gammaint_solver,
        gammaint_preconditioner=gammaint_preconditioner,
        gammaint_smoother=gammaint_smoother,
        gammaint_tolerance=gammaint_tolerance,
        gammaint_rel_tol=gammaint_rel_tol,
        relaxation_reThetat=relaxation_reThetat,
        relaxation_gammaInt=relaxation_gammaInt,
    )
    with open(output_path, "w") as file:
        file.write(content)


def fv_schemes_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the fvSchemes file for OpenFOAM simulation based on the provided
    settings.

    Args:
        setup (Settings): The simulation settings json file.
        output_path (Path): The path to save the fvSchemes file.
    """
    fv_schemes_setup = setup.simulation_settings.get("FvSchemes", {})

    time_schemes = fv_schemes_setup.get("ddtSchemes", {})
    time_scheme_default = time_schemes.get("Default", "steadyState")

    div_schemes = fv_schemes_setup.get("DivSchemes", {})
    div_scheme_default = div_schemes.get("Default", "none")
    div_scheme_U = div_schemes.get("DivPhiU", "Gauss upwind")
    div_scheme_phi_k = div_schemes.get("DivPhiK", "Gauss upwind")
    div_scheme_phi_epsilon = div_schemes.get("DivPhiEpsilon", "Gauss upwind")
    div_scheme_phi_omega = div_schemes.get("DivPhiOmega", "Gauss upwind")
    div_scheme_phi_nuTilda = div_schemes.get("DivPhiNuTilda", "Gauss upwind")

    grad_schemes = fv_schemes_setup.get("GradSchemes", {})
    grad_scheme_default = grad_schemes.get("Default", "Gauss linear")

    laplacian_schemes = fv_schemes_setup.get("LaplacianSchemes", {})
    laplacian_scheme_default = laplacian_schemes.get(
        "Default", "Gauss linear corrected")

    interpolation_schemes = fv_schemes_setup.get("InterpolationSchemes", {})
    interpolation_scheme_default = interpolation_schemes.get("Default", "linear")

    sn_grad_schemes = fv_schemes_setup.get("SnGradSchemes", {})
    sn_grad_scheme_default = sn_grad_schemes.get("Default", "corrected")

    wall_dist_schemes = fv_schemes_setup.get("WallDistSchemes", {})
    wall_dist_scheme_default = wall_dist_schemes.get("Default", "meshWave")

    content = generate_fv_schemes_dict(
        time_scheme=time_scheme_default,
        grad_scheme=grad_scheme_default,
        div_scheme_default=div_scheme_default,
        div_scheme_U=div_scheme_U,
        div_scheme_phi_k=div_scheme_phi_k,
        div_scheme_phi_epsilon=div_scheme_phi_epsilon,
        div_scheme_phi_omega=div_scheme_phi_omega,
        div_scheme_phi_nuTilda=div_scheme_phi_nuTilda,
        laplacian_scheme=laplacian_scheme_default,
        interpolation_scheme=interpolation_scheme_default,
        sn_grad_scheme=sn_grad_scheme_default,
        wall_dist_scheme=wall_dist_scheme_default
    )
    with open(output_path, "w") as file:
        file.write(content)


def decompose_par_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the decomposeParDict file for OpenFOAM simulation based on the
    provided settings.

    Args:
        setup (Settings): The simulation settings json file.
        output_path (Path): The path to save the decomposeParDict file.
    """
    number_of_subdomains = setup.simulation_settings.get(
        "Decomposition", {}).get("NumberOfSubdomains", 4)
    method = setup.simulation_settings.get("Decomposition", {}).get(
        "Method", "scotch")

    content = generate_decompose_par_dict(
        number_of_subdomains=number_of_subdomains,
        method=method
    )
    with open(output_path, "w") as file:
        file.write(content)


def surface_feature_extract_dict(setup: Settings, output_path: Path) -> None:
    """
    Generate surfaceFeatureExtractDict for extracting sharp edges (LE/TE).

    Args:
        setup (Settings): The simulation settings json file.
        output_path (Path): The path to save the surfaceFeatureExtractDict file.
    """
    snappy = setup.mesh_settings.get("SnappyHexMesh", {})
    included_angle = snappy.get("FeatureIncludedAngle", 150)

    content = generate_surface_feature_extract_dict(
        included_angle=included_angle
    )
    with open(output_path, "w") as f:
        f.write(content)


def extrude_mesh_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the extrudeMeshDict file for OpenFOAM simulation based on the provided
    settings. Thickness is computed from the bounding box z-span.

    Args:
        setup (Settings): The simulation settings json file.
        output_path (Path): The path to save the extrudeMeshDict file.
    """
    extrude_setup = setup.simulation_settings.get("ExtrudeMesh", {})
    bbox = setup.mesh_settings.get("BoundingBox", {})

    # Calculate thickness from bounding box
    z_min = float(bbox.get("ZMin", -0.5))
    z_max = float(bbox.get("ZMax", 0.5))
    thickness = abs(z_max - z_min)

    construct_from = extrude_setup.get("ConstructFrom", "patch")
    source_case = extrude_setup.get("SourceCase", ".")
    source_patches = extrude_setup.get("SourcePatches", ["front"])
    flip_normals = extrude_setup.get("FlipNormals", False)
    exposed_patch_name = extrude_setup.get("ExposedPatchName", "back")

    extrude_model = extrude_setup.get("ExtrudeModel", "linearNormal")
    n_layers = extrude_setup.get("nLayers", 1)
    expansion_ratio = extrude_setup.get("ExpansionRatio", 1.0)

    preserve_patches = extrude_setup.get("PreservePatches", False)
    merge_patch_faces = extrude_setup.get("MergePatchFaces", False)
    merge_faces = extrude_setup.get("MergeFaces", False)

    content = generate_extrude_mesh_dict(
        construct_from=construct_from,
        source_case=source_case,
        source_patches=source_patches,
        flip_normals=flip_normals,
        exposed_patch_name=exposed_patch_name,
        extrude_model=extrude_model,
        n_layers=n_layers,
        expansion_ratio=expansion_ratio,
        thickness=thickness,
        preserve_patches=preserve_patches,
        merge_patch_faces=merge_patch_faces,
        merge_faces=merge_faces,
    )
    with open(output_path, "w") as file:
        file.write(content)


def add_force_coeffs_dict_from_bc(
    control_dict_path: Path,
    bc: BoundaryConditions,
    chord: float,
    airfoil_patch: str = "airfoil",
    span: float = 1.0,
    cofr_x_frac: float = 0.25,
    angle_of_attack_deg: float = 0.0
) -> None:
    """
    Append a minimal forceCoeffs function object to an existing controlDict
    file, using values from a BoundaryConditions object.
    CofR is set at 0.25 chord, rotated by angle of attack.

    Args:
        control_dict_path (Path): Path to controlDict file.
        bc (BoundaryConditions): Boundary conditions object with velocity and density.
        chord (float): Airfoil chord length.
        airfoil_patch (str): Name of airfoil patch.
        span (float): Span length. Default: 1.0
        cofr_x_frac (float): Fraction of chord for center of rotation x-coord.
                            Default: 0.25 (quarter-chord)
        angle_of_attack_deg (float): Angle of attack in degrees. Default: 0.0
    """
    velocity = bc.velocity
    if hasattr(velocity, "__len__"):
        mag_u_inf = float(np.linalg.norm(velocity))
    else:
        mag_u_inf = float(velocity)
    a_ref = chord * span
    density = float(bc.density)

    # Calculate CofR at 0.25 chord, rotated by AoA
    alpha_rad = np.radians(angle_of_attack_deg)
    cofr_x = chord * cofr_x_frac * np.cos(alpha_rad)
    cofr_y = chord * cofr_x_frac * np.sin(alpha_rad)
    cofr = [cofr_x, cofr_y, 0]

    lift_dir = [0, 1, 0]
    drag_dir = [1, 0, 0]
    pitch_axis = [0, 0, 1]

    with open(control_dict_path, "a") as f:
        f.write(
            f"""
functions
{{
    forceCoeffs
    {{
        type            forceCoeffs;
        libs            ("libforces.so");
        patches         ({airfoil_patch});
        rho         rhoInf;
        rhoInf          {density};
        log             yes;
        writeControl    timeStep;
        writeInterval   5;
        CofR            ({' '.join(str(x) for x in cofr)});
        liftDir         ({' '.join(str(x) for x in lift_dir)});
        dragDir         ({' '.join(str(x) for x in drag_dir)});
        pitchAxis       ({' '.join(str(x) for x in pitch_axis)});
        magUInf         {mag_u_inf};
        lRef            {chord};
        Aref            {a_ref};
        writeFields     no;
    }}
}}
"""
        )
