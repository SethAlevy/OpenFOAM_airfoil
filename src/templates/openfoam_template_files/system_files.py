"""
OpenFOAM system directory file templates.

This module contains functions to generate OpenFOAM system directory files:
- controlDict
- fvSolution
- fvSchemes
- decomposeParDict
- extrudeMeshDict
- surfaceFeatureExtractDict
"""


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


def build_solver_section(
        solver: str,
        preconditioner: str,
        smoother: str,
        tolerance: float,
        rel_tol: float
) -> str:
    """
    Build a solver section with appropriate smoother/preconditioner.
    
    Args:
        solver (str): Solver type.
        preconditioner (str): Preconditioner type.
        smoother (str): Smoother type.
        tolerance (float): Solver tolerance.
        rel_tol (float): Relative tolerance.
    """
    section = f"    {{\n        solver          {solver};\n"
    if preconditioner is not None:
        section += f"        preconditioner  {preconditioner};\n"
    if smoother is not None:
        section += f"        smoother        {smoother};\n"
    section += f"        tolerance       {tolerance};\n"
    section += f"        relTol          {rel_tol};\n    }}"
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
    Generate fvSolution file content suitable for an airfoil simulation with
    essential parameters and solver sections.

    Supports both 'smoother' and 'preconditioner' keywords - only includes the
    one that is provided.

    Args:
        p_solver (str): Pressure solver type.
        p_preconditioner (str): Pressure preconditioner type.
        p_smoother (str): Pressure smoother type.
        p_tolerance (float): Pressure solver tolerance.
        p_rel_tol (float): Pressure solver relative tolerance.
        U_solver (str): Velocity solver type.
        U_preconditioner (str): Velocity preconditioner type.
        U_smoother (str): Velocity smoother type.
        U_tolerance (float): Velocity solver tolerance.
        U_rel_tol (float): Velocity solver relative tolerance.
        turb_solver (str): Turbulence parameters (k, omega, epsilon, nuTilda) solver type.
        turb_preconditioner (str): Turbulence parameters (k, omega, epsilon, nuTilda) 
            preconditioner type.
        turb_smoother (str): Turbulence parameters (k, omega, epsilon, nuTilda) 
            smoother type.
        turb_tolerance (float): Turbulence parameters (k, omega, epsilon, nuTilda) solver
            tolerance.
        turb_rel_tol (float): Turbulence parameters (k, omega, epsilon, nuTilda) solver
            relative tolerance.
        rethetat_solver (str): Rethetat solver type.
        rethetat_preconditioner (str): Rethetat preconditioner type.
        rethetat_smoother (str): Rethetat smoother type.
        rethetat_tolerance (float): Rethetat solver tolerance.
        rethetat_rel_tol (float): Rethetat solver relative tolerance.
        gammaint_solver (str): GammaInt solver type.
        gammaint_preconditioner (str): GammaInt preconditioner type.
        gammaint_smoother (str): GammaInt smoother type.
        gammaint_tolerance (float): GammaInt solver tolerance.
        gammaint_rel_tol (float): GammaInt solver relative tolerance.
        n_non_ortho_correctors (int): Number of non-orthogonal correctors for SIMPLE.
        relaxation_p (float): Relaxation factor for pressure.
        relaxation_U (float): Relaxation factor for velocity.
        relaxation_k (float): Relaxation factor for turbulent kinetic energy.
        relaxation_omega (float): Relaxation factor for specific dissipation rate.
        relaxation_epsilon (float): Relaxation factor for dissipation rate.
        relaxation_nuTilda (float): Relaxation factor for turbulent viscosity.
        relaxation_reThetat (float): Relaxation factor for ReThetat.
        relaxation_gammaInt (float): Relaxation factor for gammaInt.
        cache (str): Cache settings.

    Returns:
        str: The filled fvSolution content.
    """

    p_section = build_solver_section(
        p_solver,
        p_preconditioner,
        p_smoother,
        p_tolerance,
        p_rel_tol
    )
    U_section = build_solver_section(
        U_solver,
        U_preconditioner,
        U_smoother,
        U_tolerance,
        U_rel_tol
    )
    turb_section = build_solver_section(
        turb_solver,
        turb_preconditioner,
        turb_smoother,
        turb_tolerance,
        turb_rel_tol
    )
    rethetat_section = build_solver_section(
        rethetat_solver,
        rethetat_preconditioner,
        rethetat_smoother,
        rethetat_tolerance,
        rethetat_rel_tol
    )
    gammaint_section = build_solver_section(
        gammaint_solver,
        gammaint_preconditioner,
        gammaint_smoother,
        gammaint_tolerance,
        gammaint_rel_tol
    )

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
        div_scheme_default (str): Default divergence scheme (fall back if not specified).
        div_scheme_U (str): Divergence scheme for velocity.
        div_scheme_phi_k (str): Divergence scheme for turbulent kinetic energy.
        div_scheme_phi_epsilon (str): Divergence scheme for dissipation rate.
        div_scheme_phi_omega (str): Divergence scheme for specific dissipation
                                   rate.
        div_scheme_phi_nuTilda (str): Divergence scheme for turbulent viscosity.
        div_scheme_phi_gammaInt (str): Divergence scheme for gammaInt.
        div_scheme_phi_reThetat (str): Divergence scheme for ReThetat.
        div_scheme_nu_eff (str): Divergence scheme for effective viscosity.
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


def generate_surface_feature_extract_dict(included_angle: float = 150) -> str:
    """
    Generate surfaceFeatureExtractDict for extracting sharp edges (LE/TE).

    Args:
        included_angle (float): Angle threshold for edge extraction.

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

    Args:
        construct_from (str): Construction method.
        source_case (str): Source case path.
        source_patches (list[str]): Source patch names.
        flip_normals (bool): Flip normal direction.
        exposed_patch_name (str): Name for exposed patch.
        extrude_model (str): Extrusion model type.
        n_layers (int): Number of extrusion layers.
        expansion_ratio (float): Layer expansion ratio.
        thickness (float): Total extrusion thickness.
        preserve_patches (bool): Preserve patches.
        merge_patch_faces (bool): Merge patch faces.
        merge_faces (bool): Merge internal faces.

    Returns:
        str: The filled extrudeMeshDict content.
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
