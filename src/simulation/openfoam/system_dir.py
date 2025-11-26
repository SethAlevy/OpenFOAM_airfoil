from pathlib import Path
from templates.initial_settings_template import Settings


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
        run_time_modifiable: str = "true"
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
    p_solver = fv_solution_setup.get("p_solver", "PCG")
    p_preconditioner = fv_solution_setup.get("p_preconditioner", "DIC")
    p_tolerance = fv_solution_setup.get("p_tolerance", 1e-6)
    p_rel_tol = fv_solution_setup.get("p_rel_tol", 0)
    u_solver = fv_solution_setup.get("U_solver", "smoothSolver")
    u_smoother = fv_solution_setup.get("U_smoother", "symGaussSeidel")
    u_tolerance = fv_solution_setup.get("U_tolerance", 1e-5)
    u_rel_tol = fv_solution_setup.get("U_rel_tol", 0)
    turb_solver = fv_solution_setup.get("turb_solver", "smoothSolver")
    turb_smoother = fv_solution_setup.get("turb_smoother", "symGaussSeidel")
    turb_tolerance = fv_solution_setup.get("turb_tolerance", 1e-5)
    turb_rel_tol = fv_solution_setup.get("turb_rel_tol", 0)
    n_non_ortho_correctors = fv_solution_setup.get("n_non_ortho_correctors", 2)
    relaxation_p = fv_solution_setup.get("relaxation_p", 0.3)
    relaxation_U = fv_solution_setup.get("relaxation_U", 0.7)
    relaxation_k = fv_solution_setup.get("relaxation_k", 0.7)
    relaxation_omega = fv_solution_setup.get("relaxation_omega", 0.7)
    relaxation_epsilon = fv_solution_setup.get("relaxation_epsilon", 0.7)
    relaxation_nuTilda = fv_solution_setup.get("relaxation_nuTilda", 0.7)
    cache = fv_solution_setup.get("cache", "grad(U);")

    content = generate_fv_solution_dict(
        p_solver=p_solver,
        p_preconditioner=p_preconditioner,
        p_tolerance=p_tolerance,
        p_rel_tol=p_rel_tol,
        U_solver=u_solver,
        U_smoother=u_smoother,
        U_tolerance=u_tolerance,
        U_rel_tol=u_rel_tol,
        turb_solver=turb_solver,
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
        cache=cache
    )
    with open(output_path, "w") as file:
        file.write(content)


def generate_fv_solution_dict(
        p_solver: str = "PCG",
        p_preconditioner: str = "DIC",
        p_tolerance: float = 1e-6,
        p_rel_tol: float = 0,
        U_solver: str = "smoothSolver",
        U_smoother: str = "symGaussSeidel",
        U_tolerance: float = 1e-5,
        U_rel_tol: float = 0,
        turb_solver: str = "smoothSolver",
        turb_smoother: str = "symGaussSeidel",
        turb_tolerance: float = 1e-5,
        turb_rel_tol: float = 0,
        n_non_ortho_correctors: int = 2,
        relaxation_p: float = 0.3,
        relaxation_U: float = 0.7,
        relaxation_k: float = 0.7,
        relaxation_omega: float = 0.7,
        relaxation_epsilon: float = 0.7,
        relaxation_nuTilda: float = 0.7,
        cache: str = "grad(U);"
) -> str:
    """
    Generate fvSolution file content suitable for an airfoil simulation with essential
    parameters.
    
    Args:
        p_solver (str): Solver for pressure.
        p_preconditioner (str): Preconditioner for pressure solver.
        p_tolerance (float): Tolerance for pressure solver.
        p_rel_tol (float): Relative tolerance for pressure solver.
        U_solver (str): Solver for velocity.
        U_smoother (str): Smoother for velocity solver.
        U_tolerance (float): Tolerance for velocity solver.
        U_rel_tol (float): Relative tolerance for velocity solver.
        turb_solver (str): Solver for turbulence quantities.
        turb_smoother (str): Smoother for turbulence solver.
        turb_tolerance (float): Tolerance for turbulence solver.
        turb_rel_tol (float): Relative tolerance for turbulence solver.
        n_non_ortho_correctors (int): Number of non-orthogonal correctors.
        relaxation_p (float): Relaxation factor for pressure.
        relaxation_U (float): Relaxation factor for velocity.
        relaxation_k (float): Relaxation factor for turbulent kinetic energy.
        relaxation_omega (float): Relaxation factor for specific dissipation rate.
        relaxation_epsilon (float): Relaxation factor for dissipation rate.
        relaxation_nuTilda (float): Relaxation factor for turbulent viscosity.
        cache (str): Cache settings.

    Returns:
        str: The filled fvSolution content.
    """
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
    {{
        solver          {p_solver};
        preconditioner  {p_preconditioner};
        tolerance       {p_tolerance};
        relTol          {p_rel_tol};
    }}
    U
    {{
        solver          {U_solver};
        smoother        {U_smoother};
        tolerance       {U_tolerance};
        relTol          {U_rel_tol};
    }}
    "(k|omega|epsilon|nuTilda)"
    {{
        solver          {turb_solver};
        smoother        {turb_smoother};
        tolerance       {turb_tolerance};
        relTol          {turb_rel_tol};
    }}
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
    fv_schemes_setup = setup.simulation_settings.get("Schemes", {})
    time_scheme = fv_schemes_setup.get("TimeScheme", "steadyState")
    grad_scheme = fv_schemes_setup.get("GradScheme", "Gauss linear")
    div_scheme_U = fv_schemes_setup.get("DivScheme_U", "Gauss upwind")
    div_scheme_phi_k = fv_schemes_setup.get("DivScheme_phi_k", "Gauss upwind")
    div_scheme_phi_epsilon = fv_schemes_setup.get("DivScheme_phi_epsilon", "Gauss upwind")
    laplacian_scheme = fv_schemes_setup.get("LaplacianScheme", "Gauss linear corrected")
    interpolation_scheme = fv_schemes_setup.get("InterpolationScheme", "linear")
    sn_grad_scheme = fv_schemes_setup.get("SnGradScheme", "corrected")

    content = generate_fv_schemes_dict(
        time_scheme=time_scheme,
        grad_scheme=grad_scheme,
        div_scheme_U=div_scheme_U,
        div_scheme_phi_k=div_scheme_phi_k,
        div_scheme_phi_epsilon=div_scheme_phi_epsilon,
        laplacian_scheme=laplacian_scheme,
        interpolation_scheme=interpolation_scheme,
        sn_grad_scheme=sn_grad_scheme
    )
    with open(output_path, "w") as file:
        file.write(content)


def generate_fv_schemes_dict(
        time_scheme: str = "steadyState",
        grad_scheme: str = "Gauss linear",
        div_scheme_U: str = "Gauss upwind",
        div_scheme_phi_k: str = "Gauss upwind",
        div_scheme_phi_epsilon: str = "Gauss upwind",
        laplacian_scheme: str = "Gauss linear corrected",
        interpolation_scheme: str = "linear",
        sn_grad_scheme: str = "corrected"
) -> str:
    """
    Generate fvSchemes file content with essential parameters.

    Args:
        time_scheme (str): Time discretization scheme.
        grad_scheme (str): Gradient scheme.
        div_scheme_U (str): Divergence scheme for velocity.
        div_scheme_phi_k (str): Divergence scheme for turbulent kinetic energy.
        div_scheme_phi_epsilon (str): Divergence scheme for dissipation rate.
        laplacian_scheme (str): Laplacian scheme.
        interpolation_scheme (str): Interpolation scheme.
        sn_grad_scheme (str): SnGrad scheme.

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
    div(phi,U)      {div_scheme_U};
    div(phi,k)      {div_scheme_phi_k};
    div(phi,epsilon) {div_scheme_phi_epsilon};
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


def surface_feature_extract_dict(setup: Settings, output_path: Path) -> None:
    """
    Fill the surfaceFeatureExtractDict file for OpenFOAM simulation based on the provided

    Args:
        setup (Settings): The simulation settings.
        output_path (Path): The path to save the surfaceFeatureExtractDict file.
    """
    sfe_settings = setup.simulation_settings.get("SurfaceFeatureExtract", {})
    stl_file = sfe_settings.get("StlFile", "airfoil.stl")
    extraction_method = sfe_settings.get("ExtractionMethod", "extractFromSurface")
    included_angle = sfe_settings.get("IncludedAngle", 150)
    write_obj = sfe_settings.get("WriteObj", "yes")

    content = generate_surface_feature_extract_dict(
        stl_file=stl_file,
        extraction_method=extraction_method,
        included_angle=included_angle,
        write_obj=write_obj
    )
    with open(output_path, "w") as file:
        file.write(content)


def generate_surface_feature_extract_dict(
        stl_file: str = "airfoil.stl",
        extraction_method: str = "extractFromSurface",
        included_angle: int = 150,
        write_obj: str = "yes"
) -> str:
    """
    Generate surfaceFeatureExtractDict file content.
    
    Args:
        stl_file (str): The STL file name.
        extraction_method (str): The extraction method.
        included_angle (int): The included angle.
        write_obj (str): Whether to write OBJ file.

    Returns:
        str: The filled surfaceFeatureExtractDict content.
    """
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      surfaceFeatureExtractDict;
}}

{stl_file}
{{
    extractionMethod    {extraction_method};
    extractFromSurfaceCoeffs
    {{
        includedAngle   {included_angle};
    }}
    writeObj            {write_obj};
}}
"""
