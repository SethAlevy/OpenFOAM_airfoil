from pathlib import Path
from templates.initial_settings_template import Settings
from templates.boundary_conditions.boundary_condition_template import BoundaryCondition
import numpy as np


def control_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the controlDict file for OpenFOAM simulation based on the provided settings.
    Then write the file.

    Args:
        setup (Settings): The simulation settings.
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


def generate_control_dict(
        solver: str = "simpleFoam",
        start_from: str = "startTime",
        start_time: float = 0.0,
        stop_at: str = "endTime",
        end_time: float = 1000.0,
        delta_t: float = 1.0,
        write_control: str = "timeStep",
        write_interval: int = 100,
        purge_write: int = 0,
        write_format: str = "ascii",
        write_precision: int = 6,
        write_compression: str = "off",
        time_format: str = "general",
        time_precision: int = 6,
        run_time_modifiable: str = "true",
) -> str:
    """
    Generate controlDict file content with essential parameters.

    Args:
        solver (str): The solver to use.
        start_from (str): Time to start from.
        start_time (float): The start time.
        stop_at (str): Condition to stop at.
        end_time (float): The end time.
        delta_t (float): Time step size.
        write_control (str): Write control type.
        write_interval (int): Interval for writing data.
        purge_write (int): Number of time directories to keep.
        write_format (str): Format for writing data.
        write_precision (int): Precision for writing data.
        write_compression (str): Compression for writing data.
        time_format (str): Format for time representation.
        time_precision (int): Precision for time representation.
        run_time_modifiable (str): Whether runtime modifiable is enabled.

    Returns:
        str: The filled controlDict content.

    """
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}}

application     {solver};

startFrom       {start_from};
startTime       {start_time};
stopAt          {stop_at};
endTime         {end_time};
deltaT          {delta_t};
writeControl    {write_control};
writeInterval   {write_interval};
purgeWrite      {purge_write};
writeFormat     {write_format};
writePrecision  {write_precision};
writeCompression {write_compression};
timeFormat      {time_format};
timePrecision   {time_precision};
runTimeModifiable {run_time_modifiable};
"""
    return content


def fv_solution_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the fvSolution file for OpenFOAM simulation based on the provided settings.

    Args:
        setup (Settings): The simulation settings.
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

    n_non_ortho_correctors = fv_solution_setup.get("n_non_ortho_correctors", 2)

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


def build_solver_section(solver: str, preconditioner: str, smoother: str,
                         tolerance: float, rel_tol: float) -> str:
    """Build a solver section with appropriate smoother/preconditioner."""
    section = f"    {{\n        solver          {solver};\n"
    if preconditioner is not None:
        section += f"        preconditioner  {preconditioner};\n"
    if smoother is not None:
        section += f"        smoother        {smoother};\n"
    section += f"        tolerance       {tolerance};\n        relTol          {rel_tol};\n    }}"
    return section


def generate_fv_solution_dict(
    p_solver: str = "PCG",
    p_preconditioner: str = None,
    p_smoother: str = None,
    p_tolerance: float = 1e-6,
    p_rel_tol: float = 0,
    U_solver: str = "smoothSolver",
    U_preconditioner: str = None,
    U_smoother: str = "symGaussSeidel",
    U_tolerance: float = 1e-5,
    U_rel_tol: float = 0,
    turb_solver: str = "smoothSolver",
    turb_preconditioner: str = None,
    turb_smoother: str = "symGaussSeidel",
    turb_tolerance: float = 1e-5,
    turb_rel_tol: float = 0,
    rethetat_solver: str = "smoothSolver",
    rethetat_preconditioner: str = None,
    rethetat_smoother: str = "symGaussSeidel",
    rethetat_tolerance: float = 1e-8,
    rethetat_rel_tol: float = 0.1,
    gammaint_solver: str = "smoothSolver",
    gammaint_preconditioner: str = None,
    gammaint_smoother: str = "symGaussSeidel",
    gammaint_tolerance: float = 1e-8,
    gammaint_rel_tol: float = 0.1,
    n_non_ortho_correctors: int = 2,
    relaxation_p: float = 0.3,
    relaxation_U: float = 0.7,
    relaxation_k: float = 0.7,
    relaxation_omega: float = 0.7,
    relaxation_epsilon: float = 0.7,
    relaxation_nuTilda: float = 0.7,
    relaxation_reThetat: float = 0.4,
    relaxation_gammaInt: float = 0.4,
    cache: str = "grad(U);",
) -> str:
    """
    Generate fvSolution file content suitable for an airfoil simulation with essential
    parameters, including ReThetat and gammaInt solver sections.

    Supports both 'smoother' and 'preconditioner' keywords - only includes the one
    that is provided (not None).
    """

    p_section = build_solver_section(p_solver, p_preconditioner, p_smoother,
                                     p_tolerance, p_rel_tol)
    U_section = build_solver_section(U_solver, U_preconditioner, U_smoother,
                                     U_tolerance, U_rel_tol)
    turb_section = build_solver_section(turb_solver, turb_preconditioner, turb_smoother,
                                        turb_tolerance, turb_rel_tol)
    rethetat_section = build_solver_section(rethetat_solver, rethetat_preconditioner,
                                            rethetat_smoother, rethetat_tolerance, rethetat_rel_tol)
    gammaint_section = build_solver_section(gammaint_solver, gammaint_preconditioner,
                                            gammaint_smoother, gammaint_tolerance, gammaint_rel_tol)

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}}

solvers
{{
    p
{p_section}
    U
{U_section}
    "(k|omega|epsilon|nuTilda)"
{turb_section}
    ReThetat
{rethetat_section}
    gammaInt
{gammaint_section}
}}

SIMPLE
{{
    nNonOrthogonalCorrectors {n_non_ortho_correctors};
}}

relaxationFactors
{{
    fields
    {{
        p               {relaxation_p};
    }}
    equations
    {{
        U               {relaxation_U};
        k               {relaxation_k};
        omega           {relaxation_omega};
        epsilon         {relaxation_epsilon};
        nuTilda         {relaxation_nuTilda};
        ReThetat       {relaxation_reThetat};
        gammaInt       {relaxation_gammaInt};
    }}
}}

cache
{{
    {cache}
}}"""
    return content


def fv_schemes_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the fvSchemes file for OpenFOAM simulation based on the provided settings.

    Args:
        setup (Settings): The simulation settings.
        output_path (Path): The path to save the fvSchemes file.
    """
    fv_schemes_setup = setup.simulation_settings.get("FvSchemes", {})

    time_schemes = fv_schemes_setup.get("ddtSchemes", {})
    time_scheme_default = time_schemes.get("Default", "steadyState")

    div_schemes = fv_schemes_setup.get("DivSchemes", {})
    div_scheme_default = div_schemes.get("Default", "none")
    div_scheme_U = div_schemes.get("DivPhiU", "Gauss upwind")
    div_scheme_phi_k = div_schemes.get("DivPhiK", "Gauss upwind")
    div_scheme_phi_epsilon = div_schemes.get(
        "DivPhiEpsilon", "Gauss upwind")
    div_scheme_phi_omega = div_schemes.get("DivPhiOmega", "Gauss upwind")
    div_scheme_phi_nuTilda = div_schemes.get(
        "DivPhiNuTilda", "Gauss upwind")

    grad_schemes = fv_schemes_setup.get("GradSchemes", {})
    grad_scheme_default = grad_schemes.get("Default", "Gauss linear")

    laplacian_schemes = fv_schemes_setup.get("LaplacianSchemes", {})
    laplacian_scheme_default = laplacian_schemes.get("Default", "Gauss linear corrected")

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


def generate_fv_schemes_dict(
    time_scheme: str = "steadyState",
    grad_scheme: str = "Gauss linear",
    div_scheme_default: str = "none",
    div_scheme_U: str = "Gauss upwind",
    div_scheme_phi_k: str = "Gauss upwind",
    div_scheme_phi_epsilon: str = "Gauss upwind",
    div_scheme_phi_omega: str = "Gauss upwind",
    div_scheme_phi_nuTilda: str = "Gauss upwind",
    div_scheme_phi_gammaInt: str = "Gauss upwind",
    div_scheme_phi_reThetat: str = "Gauss upwind",
    div_scheme_nu_eff: str = "Gauss linear",
    laplacian_scheme: str = "Gauss linear corrected",
    interpolation_scheme: str = "linear",
    sn_grad_scheme: str = "corrected",
    wall_dist_scheme: str = "meshWave"
) -> str:
    """
    Generate fvSchemes file content with essential parameters.

    Args:
        time_scheme (str): Time discretization scheme.
        grad_scheme (str): Gradient scheme.
        div_scheme_U (str): Divergence scheme for velocity.
        div_scheme_phi_k (str): Divergence scheme for turbulent kinetic energy.
        div_scheme_phi_epsilon (str): Divergence scheme for dissipation rate.
        div_scheme_phi_omega (str): Divergence scheme for specific dissipation rate.
        div_scheme_phi_nuTilda (str): Divergence scheme for turbulent viscosity.
        laplacian_scheme (str): Laplacian scheme.
        interpolation_scheme (str): Interpolation scheme.
        sn_grad_scheme (str): SnGrad scheme.
        wall_dist_scheme (str): Wall distance calculation method.

    Returns:
        str: The filled fvSchemes content.
    """
    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}}

ddtSchemes
{{
    default         {time_scheme};
}}

gradSchemes
{{
    default         {grad_scheme};
}}

divSchemes
{{
    default         {div_scheme_default};
    div(phi,U)      {div_scheme_U};
    div(phi,k)      {div_scheme_phi_k};
    div(phi,epsilon) {div_scheme_phi_epsilon};
    div(phi,omega)   {div_scheme_phi_omega};
    div(phi,nuTilda) {div_scheme_phi_nuTilda};
    div(phi,gammaInt) {div_scheme_phi_gammaInt};
    div(phi,ReThetat) {div_scheme_phi_reThetat};
    div((nuEff*dev2(T(grad(U))))) {div_scheme_nu_eff};
}}

laplacianSchemes
{{
    default         {laplacian_scheme};
}}

interpolationSchemes
{{
    default         {interpolation_scheme};
}}

snGradSchemes
{{
    default         {sn_grad_scheme};
}}

wallDist
{{
    method          {wall_dist_scheme};
}}
"""
    return content


def decompose_par_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the decomposeParDict file for OpenFOAM simulation based on the provided
    settings.

    Args:
        setup (Settings): The simulation settings.
        output_path (Path): The path to save the decomposeParDict file.
    """
    number_of_subdomains = setup.simulation_settings.get("Decomposition", {}).get(
        "NumberOfSubdomains", 4)
    method = setup.simulation_settings.get("Decomposition", {}).get(
        "Method", "scotch")

    content = generate_decompose_par_dict(
        number_of_subdomains=number_of_subdomains,
        method=method
    )
    with open(output_path, "w") as file:
        file.write(content)


def generate_decompose_par_dict(
        number_of_subdomains: int = 4,
        method: str = "scotch"
) -> str:
    """
    Generate a minimal decomposeParDict file content.

    Args:
        number_of_subdomains (int): Number of subdomains.
        method (str): Decomposition method.

    Returns:
        str: The filled decomposeParDict content.
    """
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}}

numberOfSubdomains {number_of_subdomains};

method          {method};

distributed     no;

roots           ();
"""


def surface_features_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the surfaceFeaturesDict file for OpenFOAM simulation based on the provided
        settings.

    Args:
        setup (Settings): The simulation settings.
        output_path (Path): The path to save the surfaceFeaturesDict file.
    """
    sfe_settings = setup.simulation_settings.get("SurfaceFeatures", {})
    stl_file = sfe_settings.get("StlFile", "airfoil.stl")
    included_angle = sfe_settings.get("IncludedAngle", 150)
    write_obj = sfe_settings.get("WriteObj", "yes")

    content = generate_surface_features_dict(
        stl_file=stl_file,
        included_angle=included_angle,
        write_obj=write_obj
    )
    with open(output_path, "w") as file:
        file.write(content)


def surface_feature_extract_dict(setup: Settings, output_path: Path) -> None:
    """Generate surfaceFeatureExtractDict for extracting sharp edges (LE/TE)."""
    snappy = setup.mesh_settings.get("SnappyHexMesh", {})
    included_angle = snappy.get("FeatureIncludedAngle", 150)

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      surfaceFeatureExtractDict;
}}

airfoil.stl
{{
    extractionMethod    extractFromSurface;

    extractFromSurfaceCoeffs
    {{
        includedAngle   {included_angle};
    }}

    subsetFeatures
    {{
        nonManifoldEdges    no;
        openEdges           yes;
    }}

    trimFeatures
    {{
        minElem 5;
    }}

    writeObj            yes;
    writeVTK            no;
}}
"""
    with open(output_path, "w") as f:
        f.write(content)


def add_force_coeffs_dict_from_bc(
    control_dict_path: Path,
    bc: BoundaryCondition,
    chord: float,
    airfoil_patch: str = "airfoil",
    span: float = 1.0,
    cofr_x_frac: float = 0.25,
    angle_of_attack_deg: float = 0.0
) -> None:
    """
    Append a minimal forceCoeffs function object to an existing controlDict file,
    using values from a BoundaryConditions object.
    CofR is set at 0.25 chord, rotated by angle of attack.
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


def extrude_mesh_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the extrudeMeshDict file for OpenFOAM simulation based on the provided settings.
    Thickness is computed from the bounding box z-span.

    Args:
        setup (Settings): The simulation settings.
        output_path (Path): The path to save the extrudeMeshDict file.
    """
    extrude_setup = setup.simulation_settings.get("ExtrudeMesh", {})
    bbox = setup.mesh_settings.get("BoundingBox", {})

    # Calculate thickness from bounding box
    z_min = float(bbox.get("ZMin", -0.5))
    z_max = float(bbox.get("ZMax", 0.5))
    thickness = abs(z_max - z_min)

    # Required/general entries
    construct_from = extrude_setup.get("ConstructFrom", "patch")
    source_case = extrude_setup.get("SourceCase", ".")
    source_patches = extrude_setup.get("SourcePatches", ["front"])
    flip_normals = extrude_setup.get("FlipNormals", False)
    exposed_patch_name = extrude_setup.get("ExposedPatchName", "back")

    # Extrusion controls
    extrude_model = extrude_setup.get("ExtrudeModel", "linearNormal")
    n_layers = extrude_setup.get("nLayers", 1)
    expansion_ratio = extrude_setup.get("ExpansionRatio", 1.0)
    # Override thickness from setup if explicitly provided
    thickness = extrude_setup.get("Thickness", thickness)

    # Patch/face handling
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


def generate_extrude_mesh_dict(
    construct_from: str = "patch",
    source_case: str = ".",
    source_patches: list[str] | None = None,
    flip_normals: bool = False,
    exposed_patch_name: str = "back",
    extrude_model: str = "linearNormal",
    n_layers: int = 1,
    expansion_ratio: float = 1.0,
    thickness: float = 1.0,
    preserve_patches: bool = False,
    merge_patch_faces: bool = False,
    merge_faces: bool = False,
) -> str:
    """
    Generate extrudeMeshDict file content (OpenFOAM 2406-compatible).
    """
    if source_patches is None:
        source_patches = ["front"]

    patches_str = " ".join(f'"{p}"' for p in source_patches)

    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      extrudeMeshDict;
}}

constructFrom      {construct_from};
sourceCase         "{source_case}";
sourcePatches      ({patches_str});
flipNormals        {"true" if flip_normals else "false"};
exposedPatchName   {exposed_patch_name};

extrudeModel       {extrude_model};
nLayers            {n_layers};
expansionRatio     {expansion_ratio};

linearNormalCoeffs
{{
    thickness       {thickness};
}}

preservePatches    {"true" if preserve_patches else "false"};
mergePatchFaces    {"true" if merge_patch_faces else "false"};
mergeFaces         {"true" if merge_faces else "false"};
"""
