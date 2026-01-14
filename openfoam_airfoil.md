# OpenFOAM Airfoil Simulator

This repository is dedicated to performing 2D airfoil simulations using OpenFOAM (currently OpenFOAM v2406) and post-processing the results. There are multiple functions to generate, load, and export airfoils, create meshes, prepare boundary conditions, and set up the complete case structure. Predefined environments and scripts allow you to launch full analyses for several variants with post-processing and comparisons.

# Getting Started

To avoid setup difficulties and ensure repeatability of the analyses, this repository provides predefined environments with only a few prerequisites. The Python environment is defined using Poetry, while for OpenFOAM a Docker image has been created.

The only strict requirement is Docker, which can be downloaded and installed from the [official page](https://www.docker.com/get-started/). Once the installation and configuration process is finished, the image may be built. The commands below are shown for PowerShell on Windows, but are directly transferable to Linux and macOS terminals.

While inside the repository, launch Docker Desktop and execute the following command in the terminal:

```powershell
docker build --progress=plain -t openfoam-airfoil .
```

This builds the image specified inside the [`Dockerfile`](https://github.com/SethAlevy/OpenFOAM_airfoil/blob/main/Dockerfile). The image is based on Ubuntu 22.04 and includes all system dependencies, Python 3.11 with Poetry, and OpenFOAM v2406 with cfMesh. The `--progress=plain` parameter is optional and provides additional logs to track potential build errors. All scripts inside the repository may be executed and used within it. Usually it is useful to launch it with an attached local directory where the OpenFOAM cases and JSON files with simulation settings will be located. This allows access to the files from the file explorer directly. This can be done with the following command:

```powershell
docker run -it -v "/input/dir/path:/app/case_dir" openfoam-airfoil
```

Now the repository and the case directory may be found under the `/app` path. For example, the main case running script will be under the path: `/app/src/simulation/run/run_case.sh`. You may launch it using one of the example inputs with the following commands:

```powershell
cd /app/src/simulation/run/
bash ./run_case.sh --working-dir /app/case_dir --setup-file /app/examples/naca4415_study/aoa_5.json --case-name aoa_5
```

This example launches a full workflow for a generated NACA4415 airfoil with an angle of attack of 5° and Reynolds number of 1,000,000 set in the boundary conditions.

The scripts are prepared to utilize Python scripts via the Poetry environment and it can also be used manually inside the container, but sometimes it is useful to set it up locally. A concise installation guide is available in its [documentation](https://python-poetry.org/docs/). One of the key project files for Poetry is [`pyproject.toml`](https://github.com/SethAlevy/OpenFOAM_airfoil/blob/main/pyproject.toml) in the repository root. It defines required libraries and their versions. To ensure compatibility and proper operation, versions are pinned.

To build the environment, execute:

```powershell
poetry install
```

Now it may be used through the code editor for notebooks and scripts or through the terminal:

```powershell
poetry run python src/simulation/preparation/prepare_case.py --working-dir . --setup-file examples/naca4415_study/aoa_5.json --case-name case_name
```

Lastly, you may desire to visualize the simulation results or any partial files. For this purpose it is useful to have ParaView installed. It is an open-source tool recommended for post-processing and visualization of OpenFOAM results. It may be downloaded from the [official site](https://www.paraview.org/download/).

# Airfoil

In general, an airfoil is an aerodynamic shape specially designed to generate lift (due to pressure difference on its lower and upper surfaces) and minimum drag while moving through the air. They can be found in wings or propellers as their cross-sections.

## Airfoil Terminology

There is a set of common parameters that are used to describe and understand the geometry and aerodynamic behavior of most airfoils. The most important terms are listed below:

 - **Leading edge** - the foremost point of the airfoil that encounters the incoming flow.
 - **Trailing edge** - the rearmost point of the airfoil where the upper and lower surfaces meet and the flow leaves the airfoil.
 - **Chord** - a straight line joining the leading edge and trailing edge of the airfoil. The chord length is one of the most important geometric parameters, commonly used as a reference dimension in aerodynamic equations.
 - **Camber** - a curve equivalent to the geometrical centerline, located midway between the upper and lower surfaces, describing the airfoil's curvature.
 - **Thickness** - the distance between the upper and lower surfaces. Thickness is most commonly measured perpendicular to the camber line, though it may also be defined perpendicular to the chord line. The maximum thickness is a frequently used geometric parameter.
 - **Upper and lower surfaces** - create the airfoil's outer geometry connecting the leading and trailing edges. Typically the upper surface encounters lower static pressure than the lower, and sometimes they are referred to as suction and pressure surfaces. The pressure difference between them generates lift.
 - **Aerodynamic center** - a point along the chord where the pitching moment stays unchanged regardless of the change in angle of attack and fluid speed.
 - **Center of pressure** - a point where the average pressure force is considered to act. Its location varies with changing conditions.
 - **Angle of attack** - the angle between the relative wind vector and the chord line.
 - **Lift** - the component of the aerodynamic force acting perpendicular to the incoming flow direction.
 - **Drag** - the component of the aerodynamic force acting parallel to the incoming flow direction.
 
Place for plot

## NACA Airfoils 

The National Advisory Committee for Aeronautics (NACA) developed and tested a series of airfoils in the first half of the 20th century. They were designed to provide systematic, well-defined geometries that could be easily reproduced, analyzed, and compared in experiments and calculations. Their most important feature is the description of key geometrical parameters through mathematical equations. Main parameters may be determined through the designations. Although there are modern, more capable airfoils developed, the NACA-series are widely used for education, aerodynamic studies, numerical simulations, validations, and experiments. There are two main groups of NACA airfoils.

### NACA 4-digit series

The simplest series where the geometry is described by three parameters through four digits in the form:

**NACA MPXX**

where: 
 - **M** is the maximum camber as a percentage of the chord, 
 - **P** is the position of maximum camber in tenths of the chord,
 - **XX** is the maximum thickness as a percentage of the chord. 
 
For example, NACA 4415 (4-4-15) has a maximum camber of 4% of the chord, located at 40% of the chord line from the leading edge and a maximum thickness of 15% of the chord. Another good example is the NACA 0012 (0-0-12) which has zero camber (which means the airfoil is symmetric and does not produce lift at 0° angle of attack) and a maximum thickness of 12% of the chord. The following equations allow calculation of the airfoil's geometry: 

#### Normalized parameters

Using a normalized chordwise coordinate \( x \in [0,1] \):

\[
m = \frac{M}{100}, \quad
p = \frac{P}{10}, \quad
t = \frac{XX}{100}
\]

---

#### Mean camber line

\[
y_c(x) =
\begin{cases}
\frac{m}{p^2}(2px - x^2), & x \le p \\
\frac{m}{(1-p)^2}\left[(1 - 2p) + 2px - x^2\right], & x > p
\end{cases}
\]

---

#### Thickness distribution

\[
y_t(x) =
5t\left(
0.2969\sqrt{x}
- 0.1260x
- 0.3516x^2
+ 0.2843x^3
- 0.1036x^4
\right)
\]

The coefficient \(-0.1036\) corresponds to a closed trailing edge.  
For a finite trailing-edge thickness, it may be replaced with \(-0.1015\).

---

#### Surface inclination

\[
\theta(x) = \arctan\left(\frac{dy_c}{dx}\right)
\]

\[
\frac{dy_c}{dx} =
\begin{cases}
\frac{2m}{p^2}(p - x), & x \le p \\
\frac{2m}{(1-p)^2}(p - x), & x > p
\end{cases}
\]

---

#### Upper and lower surface coordinates

\[
x_u = x - y_t \sin\theta, \quad
y_u = y_c + y_t \cos\theta
\]

\[
x_l = x + y_t \sin\theta, \quad
y_l = y_c - y_t \cos\theta
\]

---

### NACA 5-digit series

Slightly more complex in geometry and description than the 4-digit series. Their naming convention is in the form:

**NACA LPQXX** 

where: 
 - **L** after multiplying by 3/20 is the design ideal lift coefficient, 
 - **P** after multiplying by 0.05 is the location of the maximum camber in tenths of the chord, 
 - **Q** identifies the type of camber (0 if simple, 1 if reflexed), 
 - **XX** is the maximum thickness as a percentage of the chord. 
 
For example, NACA 23012 (2-3-0-12) has a design lift coefficient of 0.3, maximum camber located at 15% of the chord, maximum thickness of 12% of the chord, and a simple camber.

**Normalized parameters**

\[
c_l = 0.15L, \quad
p = 0.05P, \quad
t = \frac{XX}{100}
\]

The thickness distribution is identical to the NACA 4-digit series.

**Mean camber line (normal camber, Q = 0)**

\[
k_1 = \text{empirical constant from NACA tables}
\]

\[
y_c(x) =
\begin{cases}
\frac{k_1}{6}(x^3 - 3px^2 + p^2(3-p)x), & x \le p \\
\frac{k_1 p^3}{6}(1 - x), & x > p
\end{cases}
\]

\[
\frac{dy_c}{dx} =
\begin{cases}
\frac{k_1}{6}(3x^2 - 6px + p^2(3-p)), & x \le p \\
-\frac{k_1 p^3}{6}, & x > p
\end{cases}
\]

**Mean camber line (reflex camber, Q = 1)**

\[
k_2 = \text{empirical constant from NACA tables}
\]

\[
y_c(x) =
\begin{cases}
\frac{k_2}{6}\left(x^3 - 3px^2 + p^2(3-p)x + 0.2025(p-x)^3\right), & x \le p \\
\frac{k_2 p^3}{6}(1 - x), & x > p
\end{cases}
\]

\[
\frac{dy_c}{dx} =
\begin{cases}
\frac{k_2}{6}\left(3x^2 - 6px + p^2(3-p) - 0.6075(p-x)^2\right), & x \le p \\
-\frac{k_2 p^3}{6}, & x > p
\end{cases}
\]

**Surface coordinates**

\[
\theta = \arctan\left(\frac{dy_c}{dx}\right)
\]

\[
x_u = x - y_t \sin\theta, \quad
y_u = y_c + y_t \cos\theta
\]

\[
x_l = x + y_t \sin\theta, \quad
y_l = y_c - y_t \cos\theta
\]


### Other NACA airfoil series

Additional NACA airfoil families, such as the 6-digit and 16-digit series, were developed for more specialized aerodynamic performance, particularly with respect to pressure distribution and laminar flow control. These series use more complex definitions and are not supported by the current implementation, and therefore are not discussed further.

## UIUC Database

Maintained by the University of Illinois, the site is a widely used source of information about airfoils and other aerodynamic-related resources. At the moment it provides more than 1,600 geometries that may be downloaded in the form of a DAT file. This format allows automatic implementation into the repository knowing only the exact designation of the desired geometry.

## Repository Airfoil Tools

The repository has a set of tools that aim to make working with airfoils easier. There are three options to get a geometry: 

 - Generating a NACA 4-digit airfoil
 - Generating a NACA 5-digit airfoil
 - Importing an airfoil from the UIUC database by designation

Dedicated code may be found in `src/airfoil`, where [airfoil_base.py](https://github.com/SethAlevy/OpenFOAM_airfoil/blob/main/src/airfoil/airfoil_base.py) contains the general airfoil definitions and [airfoil.py](https://github.com/SethAlevy/OpenFOAM_airfoil/blob/main/src/airfoil/airfoil.py) contains specific classes. They may be initialized as follows:

```python
from airfoil.airfoil import NACA4, NACA5, UIUCAirfoil

naca2412 = NACA4(designation="2412", chord_length=0.6, resolution=1000)
naca23009 = NACA5(designation="23009", chord_length=1.5, resolution=1000)
b737d = UIUCAirfoil("b737d", chord_length=0.8, resolution=1000)
```

Each of these calls creates an airfoil object. Thanks to the unified code, all types of initializations may be treated in exactly the same way. Key geometry parameters are defined and may be accessed easily. These are: chord, mean camber, upper and lower surfaces, thickness, and alpha (angle of attack). Usage examples:

```python
naca2412.chord
naca2412.mean_camber_line
naca2412.upper_surface
naca2412.lower_surface
naca2412.thickness
naca2412.alpha
```

Printing out the following parameters allows you to show some important assumptions about the geometry creation in this repository. All coordinate arrays start with the point [0.0, 0.0] which corresponds to the leading edge. The airfoils will always have their leading edge at the origin of the coordinate system and the chord extends in the x-axis direction. The remaining code, mesh creation, and OpenFOAM actions will be based on this assumption. In addition, you may notice that the angle of attack is 0. An airfoil is always initialized in this position. The angle of attack is applied as a rigid-body rotation about the leading edge, without modifying the underlying geometry. It may be applied dynamically to the object as follows:

```python
b737d.set_angle_of_attack(10)
print(b737d.alpha)

b737d.set_angle_of_attack(15)
print(b737d.alpha)
```

After the airfoil has been created with the desired parameters, it is possible to visualize the current geometry with its position in space and export it as an STL file (2D in XY plane or 3D with spanwise z-axis extrusion).

```python
from pathlib import Path

b737d.plot(title="Boeing 737D Airfoil - from settings", save_path=Path("b737d.png"), show=True)
b737d.to_stl(Path("b737d_3d.stl"), dimension=3)
b737d.to_stl(Path("b737d_2d.stl"), dimension=2)
```

# Airflow Physics and Numerical Modeling

This chapter will introduce some basic concepts and definitions about airfoils and airflow that are crucial to understanding and performing aerodynamic simulations like those in the repository. Despite covering even some simple definitions, it is not intended to replace textbooks and lectures and makes use of assumptions and simplifications.

## Airflow basics

Using air as the working fluid implies certain consequences when it comes to describing the physics and mathematics. It has a set of characteristics and parameters that are used for modeling. In addition, some assumptions must be made depending on the intended application, accuracy, etc. In the simulations presented here, air is treated as:

 - A continuum fluid, meaning molecular effects are neglected and the flow is described by averaged field quantities (velocity, pressure).
 - A Newtonian fluid, where viscous stresses are linearly proportional to strain rates.
 - Incompressible, as long as the Mach number remains sufficiently low.

Under these assumptions, air is fully characterized by density (\(\rho\)) and dynamic viscosity (\(\mu\), or equivalently kinematic viscosity (\(\nu = \mu / \rho\)). Despite this, there are several important parameters and definitions that should be mentioned.

### Flow velocity and reference scales

When considering free-stream aerodynamic flows around an airfoil, there are two main factors that must be known:

 - **Free-stream velocity** - typically imposed at the inlet and far field, assumed for the given conditions.
 - **Characteristic length** - for airfoils, usually the chord length.

These reference values are not arbitrary, they define the scaling of the flow and appear explicitly in all relevant dimensionless numbers.

### Reynolds number

The Reynolds number is a dimensionless quantity used to describe flow behavior. Its value is an important parameter for viscous flow as it represents the ratio between inertial forces (instability and mixing) and viscous forces (damping of velocity gradients), and thus governs the transition between laminar (low Reynolds) and turbulent (high Reynolds) flows. It may be calculated using the following equation:

\[
\mathrm{Re} = \frac{\rho\,U\,c}{\mu} = \frac{U\,c}{\nu}
\]

where:
- ρ = fluid density
- U = free-stream velocity
- c = characteristic length (chord)
- μ = dynamic viscosity
- ν = kinematic viscosity

### Mach number

The Mach number is a dimensionless quantity representing the ratio between the local flow velocity and the speed of sound. It governs the transition between subsonic, transonic, and supersonic flows.

\[
\mathrm{Ma} = \frac{U}{a}
\]

where \(a\) is:

\[
a = \sqrt{\gamma\,R\,T}
\]

where:
- γ = ratio of specific heats (1.4 for air)
- R = specific gas constant (287.05 J/(kg·K) for air)
- T = local temperature (in K)

For truly incompressible simulations (using solvers like `simpleFoam`), the Mach number should remain below 0.3. Above this threshold, density variations (compressibility effects) begin to significantly alter the lift and drag coefficients. If your simulation requirements exceed 0.3, it is recommended to switch to a compressible solver (such as `rhoSimpleFoam`) to account for these physical changes.

### How to get fluid parameters

The two mentioned basic quantities, free-stream velocity and chord length, are usually the result of project requirements. However, to perform a numerical simulation, several other parameters of the fluid must be known (at least density, viscosity, pressure). They may be assumed, estimated based on expected environmental conditions, or taken as default values for air. Often, standard conditions for aerodynamics are air at sea level with 15 °C (288.15 K), density of 1.225 kg/m³, **dynamic viscosity of \(1.789\times10^{-5}\,\text{kg/(m·s)}\)** (kinematic viscosity \(\approx 1.46\times10^{-5}\,\text{m}^2/\text{s}\)), and pressure of 101.325 kPa (although for numerical calculations where relative pressure values matter, it is convenient to use 0 as reference).

For slightly more accurate properties at different conditions, the International Standard Atmosphere (ISA) may be used. The model provides equations to estimate temperature, pressure, and density at different altitudes assuming standard conditions at sea level.

**Air temperature at altitude:**
\[
T_h = T_0 - L\,h
\]

**Air pressure at altitude:**
\[
p_h = p_0 \left(1 - \frac{L\,h}{T_0}\right)^{\frac{g\,M}{R\,L}}
\]

**Air density at altitude:**
\[
\rho_h = \frac{p_h}{R_{\text{specific}}\,T_h}
\]

Where:
- \(L\) = Temperature lapse rate: 0.0065 [K/m]
- \(g\) = Gravity acceleration: 9.80665 [m/s²]
- \(M\) = Molar mass of air: 0.0289644 [kg/mol]
- \(R\) = Universal gas constant: 8.3144598 [J/(mol·K)]
- \(R_{\text{specific}}\) = Specific gas constant for dry air: 287.05 [J/(kg·K)]

When using conditions different from standard, you may also adjust the viscosity.

**Dynamic viscosity (Sutherland):**
\[
\mu = \mu_\text{ref} \left(\frac{T_h}{T_\text{ref}}\right)^{\tfrac32} \frac{T_\text{ref} + S}{T_h + S},
\quad \mu_\text{ref}=1.789\times 10^{-5}\ \text{kg/(m·s)},\ T_\text{ref}=288.15\,\text{K},\ S=110.4\,\text{K}
\]
**Kinematic viscosity:**
\[
\nu = \frac{\mu}{\rho_h}
\]

Where:
- \(\mu_\text{ref}\) = 1.789×10⁻⁵ [kg/(m·s)]
- \(T_\text{ref}\) = 288.15 [K]
- \(S\) = 110.4 [K]

The above equations are implemented and supported in the repository.

## Laminar, turbulent, and transitional flow

In fluid mechanics, we may distinguish two types of flow regime: laminar and turbulent. Both have their own characteristics affecting the approach taken for numerical simulations.

### Laminar flow

Laminar flow occurs when viscous forces dominate and all fluid particles move in an organized, layered form. Thus it is sometimes referred to as viscous or streamlined flow. Velocity vectors are aligned with the mean flow direction and the fluid moves in predictable, parallel paths. In numerical simulations, this allows for a relatively simple direct solution.

### Turbulent flow

Turbulent flow is characterized by chaotic variations. The particles move in paths that cross each other and in directions varying from the general velocity direction. They may form different vortex structures like swirls and eddies. The exact behavior is almost unpredictable and requires approximate methods for calculations. For typical RANS (Reynolds-Averaged Navier-Stokes) simulations, the modeling parameter of eddy viscosity is introduced, which represents the average chaotic behavior of turbulence. There are several turbulence models popular in aerodynamics that add additional equations to calculate parameters allowing description of turbulence with expected accuracy.

#### One-equation models

Simple models introducing only one additional transport equation for a turbulence quantity.

 - **Spalart-Allmaras** - a very popular and simple model dedicated to aerodynamics. It introduces a single transport equation for a modified version of the kinematic eddy viscosity. It is focused on solving the near-wall boundary layers while treating the far-field turbulence in a simplified manner. It may fail for more complex flows with strong separation or recirculation.

#### Two-equation models

More complicated models introducing two additional equations, typically solving the turbulence kinetic energy (k) and some measure of its dissipation rate or frequency.

 - **k–ε family** - a group popular in industrial applications, known for its robustness and versatility. **Turbulence kinetic energy (k)** and **dissipation rate (ε)** are solved additionally. It handles complicated turbulence well but has limitations near walls. To estimate behavior in the boundary layers, it uses wall functions instead of direct solving.

 - **k–ω family** - focused on the near-wall region, allows for higher accuracy in the boundary layer and good detection of separations. To predict turbulence, equations for **kinetic energy (k)** and **specific dissipation rate (ω)** are solved. Although it allows for high accuracy near walls, it is very sensitive to ω values in the freestream, especially for initial conditions.

 - **SST (Shear Stress Transport)** - uses the advantages of **k–ε models in free-stream regions and k–ω models near the wall**. It uses a blending function that switches between the models based on distance from wall boundaries and local flow conditions. SST is the most popular choice for aerodynamic simulations.

### Transitions

To identify the flow regime, the Reynolds number is used. For external aerodynamics on airfoils, the transition typically occurs at relatively high values, often starting around \(5 \times 10^5\), depending on surface roughness and free-stream turbulence. Considering a domain, the local fluid velocity may vary significantly, especially in the near-wall region where at the boundary the velocity approaches zero. This implies the existence of areas with laminar and areas with turbulent regimes. To obtain highly accurate results, it is crucial to correctly capture both and the transition between them. To do this, transitional turbulence models were developed.

One example is the **SST Langtry-Menter** transitional model, sometimes referred to as γ–Re(\theta\). It introduces into the SST model two additional differential equations for **intermittency (γ)** and **transition momentum-thickness Reynolds number (Re\(\theta\))**. It aims to modify turbulent transport equations to simulate laminar and laminar-to-turbulent transition behavior.

## Why airfoils generate lift

As mentioned in previous chapters, airfoils are specialized aerodynamic shapes. When airflow encounters them, the air particles are redirected, generating a total aerodynamic force. The component of this force aligned with the incoming flow direction is called drag, while the perpendicular component is called lift.

The generation of this force is a complex phenomenon that can be viewed through different physical lenses. These are not competing theories, but rather different ways of describing the same conservation laws.

### The Bernoulli perspective (conservation of energy)

The most common explanation utilizes Bernoulli’s principle, which relates pressure to velocity. Due to the airfoil's shape and angle of attack, the air is forced to move at different velocities along the upper and lower surfaces. A change in velocity results in a change in static pressure. For two points in the flow, the pressure difference can be expressed as:

\[
p_1 - p_2 = \tfrac12\,\rho\,\left(U_2^2 - U_1^2\right)
\]

### The Newtonian perspective (conservation of momentum)

Another fundamental way to explain lift is through Newton’s Third Law. As the airfoil moves through the air, it deflects the incoming flow downward (a phenomenon known as downwash). Because the airfoil exerts a force on the air to change its momentum downward, the air exerts an equal and opposite reaction force upward on the airfoil.

### Viscosity and the Kutta Condition

For a numerical simulation to produce lift, a crucial physical requirement must be satisfied: the Kutta Condition. In a theoretical, perfectly frictionless (inviscid) flow, air would try to "whip" around the sharp trailing edge from the bottom to the top, resulting in zero net lift (D'Alembert's Paradox). In reality, fluid viscosity prevents this. The air is forced to leave the sharp trailing edge smoothly, which establishes the circulation ($\Gamma$) around the airfoil and dictates the velocity differences described by Bernoulli.

### Aerodynamic coefficients

Airfoil performance depends on its shape and angle of attack. To mathematically quantify it there are three main aerodynamic coefficients in use: the drag coefficient (\(C_d\)), lift coefficient (\(C_l\)), and moment coefficient (\(C_m\)). The main scope of performing airfoil analyses, both experimental and numerical, is to derive them empirically, as they are very hard to estimate using analytical methods.

 - **Lift coefficient (\(C_l\))** - The Lift Coefficient represents the lifting capability of the airfoil relative to the dynamic pressure of the flow and the surface area. It is primarily influenced by the shape of the airfoil and the Angle of Attack (\(\alpha\)).
 - **Drag Coefficient (\(C_d\))** - The Drag Coefficient quantifies the resistance of the airfoil as it moves through the air. In CFD, this value is the sum of two components: Pressure Drag (caused by the pressure difference between the front and back) and Skin Friction Drag (caused by the air sticking to the surface due to viscosity).
 - **Moment coefficient (\(C_m\))** - The Moment Coefficient represents the aerodynamic twisting force (torque) acting on the airfoil. It is typically calculated relative to a specific point, such as the quarter-chord (25% of the chord length). This value is critical for determining the longitudinal stability of an aircraft, a negative value usually indicates a nose-down pitching moment.

## What are numerical simulations

Computational Fluid Dynamics (CFD) is a branch of fluid mechanics that uses numerical analysis and data structures to analyze and solve problems involving fluid flows. In the context of airfoil design, it serves as a "virtual wind tunnel," allowing engineers to predict how air will interact with a geometry without building physical prototypes.

While Bernoulli’s equation focuses on the conservation of energy and Newton’s laws focus on the conservation of momentum, a complete description of the flow must also satisfy the conservation of mass. In CFD, we move beyond these simplified algebraic relationships to solve the Navier-Stokes equations. These equations account for the simultaneous conservation of mass, momentum, and energy while explicitly including the effects of viscosity. This allows the simulation to capture complex phenomena like boundary layer growth and flow separation that simplified theories cannot predict.

### Navier-Stokes equations

**Continuity equation (conservation of mass)**

\[
\nabla \cdot \mathbf{u} = 0
\]

**Momentum equation**

\[
\rho \left( \frac{\partial \mathbf{u}}{\partial t} + \left( \mathbf{u} \cdot \nabla \right)\mathbf{u} \right)
= -\,\nabla p + \mu \nabla^{2}\mathbf{u} + \rho\,\mathbf{f}
\]

**Energy equation**

\[
\rho\,C_p \left( \frac{\partial T}{\partial t} + \mathbf{u} \cdot \nabla T \right)
= \nabla \cdot \left( k\,\nabla T \right) + \Phi + S_h
\]

Where:

 - \(\mathbf{u}\) – velocity vector  
 - \(p\) – static pressure  
 - \(\rho\) – density  
 - \(\mu\) – dynamic viscosity  
 - \(T\) – temperature  
 - \(C_p\) – specific heat  
 - \(k\) – thermal conductivity  
 - \(t\) – time  

These three equations are the core of CFD and describe the relationship between velocity, pressure, and density. It is worth mentioning that for low-Mach, incompressible flows the energy equation may be omitted.

### Interpretation of mathematical terms

Each element of the equations has a specific physical interpretation, thus it is worth mentioning them.

 - **Transient term (\(\partial/\partial t\))** - Represents the change in a property (velocity or temperature) over time. This is zero for steady-state simulations.
 - **Convective term (\(\mathbf{u}\cdot\nabla\))** - Describes how the flow carries a property (like momentum or heat) from one location to another.
 - **Pressure gradient (\(\nabla p\))** - The force that accelerates air from high-pressure zones toward low-pressure zones. 
 - **Diffusion/viscous terms (\(\nabla^2\))** - Describes the smoothing or spreading effect of viscosity (momentum) or thermal conductivity (heat). This creates the boundary layer near the airfoil surface.
 - **Viscous dissipation (\(\Phi\))** - The heat generated by the fluid's internal friction—critical in high-speed boundary layers.

### Approaching CFD

Although the Navier-Stokes equations are fundamental, they are not solved directly. Their formulation makes them practical only for a few very simple analytical examples. For real-world engineering applications, they are nearly impossible to solve exactly, and approximate numerical approaches must be utilized. This process involves three core elements:

 - **Spatial discretization** - The continuous 2D or 3D domain around the airfoil is divided into a grid of small elements called control volumes or cells. The shape and size of these cells are critical, the mesh must be fine enough to resolve the high velocity gradients in the boundary layer, yet coarse enough in the far-field to remain computationally efficient.

 - **Numerical method** - There are different ways to convert the calculus of the equations into a system of algebraic equations. OpenFOAM utilizes the Finite Volume Method (FVM). There are a few things specific to that method. Unlike other methods that solve equations at specific points, FVM calculates the flux (the flow of a quantity) across the faces of each cell. It is inherently conservative, it ensures that the amount of mass or momentum entering a cell is exactly equal to the amount leaving it (plus or minus any internal sources). The FVM takes the complex partial differential equations and converts them into a massive system of linear algebraic equations ($Ax=b$) that the computer can solve.

 - **Iteration and Convergence** - The solver starts with an initial guess and repeatedly calculates the flow field. With each iteration, the error (or residual) should decrease. When the error is sufficiently small, the simulation is considered converged, providing a physically consistent snapshot of the fields.

# OpenFOAM introduction

OpenFOAM is an open source software for CFD modeling. It offers a wide range of options for meshing, several solvers for different physical phenomena and built-in post-processing functions. Access to its code gives the possibility to add custom changes and functionalities for very specific applications. 

For users unfamiliar with programming, OpenFOAM may seem unintuitive and complicated. It was developed for Linux operating systems and on Windows it is typically used via WSL or Docker. Unlike other, commercial software it does not have a traditional GUI. Performing simulations with OpenFOAM requires preparing specific directory structures and text dictionary files, then launching actions by executing commands in a terminal. It offers a lot of possibilities, but has a high barrier of entry and requires the user to have some basic knowledge of what they want to achieve. This chapter covers the functionalities useful in the repository and explains how the process is structured.

## Case structure

Below is a diagram of a typical OpenFOAM case directory structure, followed by a brief explanation of each key directory and file:

```text
case_name/
├── 0/                           # Initial and boundary condition files for each field (U, p, k, omega, etc.)
│   ├── U
│   ├── p
│   └── ...
├── constant/                   # Physical properties and mesh
│   ├── polyMesh/               # Mesh files (created after meshing)
│   ├── triSurface/             # STL or FMS geometry files (for snappyHexMesh/cfMesh)
│   ├── transportProperties     # Properties of the working fluid
│   └── turbulenceProperties    # Turbulence model selection of parameters
├── system/                     # Simulation control and numerical settings
│   ├── controlDict             # Time, write, and run controls
│   ├── fvSchemes               # Discretization schemes
│   ├── fvSolution              # Solver and algorithm controls
│   ├── decomposeParDict        # Used for parallel computing
│   ├── blockMeshDict           # Base mesh settings (e.g., for snappyHexMesh, for cfMesh not required)
│   ├── meshDict                # cfMesh settings file
│   └── ...
├── postProcessing/             # (Optional) Postprocessing results
└── ...                         # Other files/scripts as needed
```

**Key elements:**
- `0/`: Contains the initial and boundary conditions for all fields. Each file inside defines the initial state of the domain in every cell for a given value. 
- `constant/`: Contains mesh (`polyMesh/` in form of points, cells, faces), geometry (`triSurface/` in the appropriate source format), and physical properties (`transportProperties`, `turbulenceProperties`).
- `system/`: Contains simulation control files (`controlDict`, `fvSchemes`, `fvSolution`, etc.).
- `postProcessing/`: (Optional) Stores results from function objects or postprocessing utilities.

This structure is required for every OpenFOAM case. The most important directories are `0/`, `constant/`, and `system/`. These must always be present for a valid case setup.

## Creating mesh

The mesh is one of the crucial elements in every simulation. It must be properly defined, with good quality elements and adjusted density in the interesting areas. As it will be the numerical representation of the domain, at this stage patches will be introduced. This chapter will focus on cfMesh which is the primary mesher of this repository. 

### Defining airfoil geometry

An airfoil geometry is complicated and hard to describe directly. For such cases, OpenFOAM supports importing it in another format. Usually 3D files in STL format are valid inputs, though cfMesh prefers composite components in FMS format. For export and conversion of airfoils there are built-in functions inside the repository (look in the Airfoil chapter). It is crucial to create the STL file in a good resolution. Depending on the chord length it may require thousands of points along the x-direction to allow the mesher to map its curvature well. 

The STL files created in the repository are in ASCII format. If you try to open one in a text editor you will see something like that:

```text
solid b'airfoil'
facet normal 0.000000 0.000000 -0.000001
  outer loop
    vertex 0.000000 0.000000 0.001000
    vertex 0.000625 0.003866 0.001000
    vertex 0.001050 0.005441 0.001000
  endloop
endfacet

...
```

Solid defines the described structure, then triangle faces are described that compose the STL.

### Semi‑2D geometries

OpenFOAM in general is a 3D oriented software. Although it has some functionalities dedicated for 2D cases it does not fully support such simulations. To overcome this a semi-2D mesh must be prepared, which may be tricky and cause unexpected errors. The approach is to create a 3D mesh with a thickness of only one cell in one of the directions, usually z-direction. 

It is very important to assume some thickness of the domain and use it consistently in all stages, from extrusion of the STL and mesher bounding box to coefficient calculations in `controlDict`.

### The domain

The volume inside which calculations are performed is the domain. It may take different forms, but for this airfoil simulations it will be a rectangular bounding box. There are no strict rules how to shape the domain besides it should be large enough to allow free-stream formation and prevent interference of the boundary conditions with airfoil effects. Here it will be usually assumed ~5 x chord inlet distance, ~10 x chord outlet distance and ~4 x chord for the top and bottom distances.

### CfMesh

cfMesh is among the most popular meshing software for OpenFOAM. It does not offer that many options like snappyHexMesh and may be problematic on very complex geometries, but is a good selection for hex-cell-based 2D mesh. To create it, here the `cartesian2DMesh` utility is used, which requires a properly prepared `meshDict` inside the `system/` directory. Here a common example is presented:

```c++
FoamFile
// header telling about the file type and location
{
    version   2.0;
    format    ascii;
    class     dictionary;
    location  "system";
    object    meshDict;
}

// definition of the domain
surfaceFile         "constant/triSurface/domain.fms";
maxCellSize         0.0399996;
locationInMesh      (-8.33 0.0);

// refinement section for objects defined inside the domain
objectRefinements
{
    airfoilBox
    {
        type            box;
        centre         (4.25 0.0 0.0);
        lengthX        9.5;
        lengthY        3.0;
        lengthZ        1;
        cellSize       0.01;
    }
    leadingEdgeSphere
    {
        type            sphere;
        centre          (0.0 0.0 0.0);
        radius          0.1;
        cellSize        0.0005;
    }

}


// refinement section for local patches inside the domain
localRefinement
{
    airfoil
    {
        cellSize 0.0015;
        refinementThickness 0.015;
    }
}

// boundary layers section for local patches
boundaryLayers
{
    patchBoundaryLayers
    {
        airfoil
        {
            nLayers              18;
            thicknessRatio       1.2;
            maxFirstLayerThickness 5e-05;
            allowedDiscontinuity 1;
        }
    }
}
```

The file consists of several dictionaries that are dedicated to adjusting the mesh in a desired manner.

Going from above we have the **definition of the domain**. In this scenario the geometry is given by the surfaceFile, next `maxCellSize` defines the maximum base-cell edge length for the whole meshing domain and locationInMesh indicates which part of the geometry should be filled with mesh.

```text
// Path to the FMS where airfoil and boundary areas in XY directions are provided.
surfaceFile         "constant/triSurface/domain.fms"; 
// This will be the cell-edge size from which all further refinements will begin. It is important that the value should be slightly uneven to avoid cfMesh division errors when it is the exact divider of the domain edge.
maxCellSize         0.0399996; 
// Point in 2D space indicating the area to generate mesh inside, as it is semi-2D x, y, and z coordinate must be provided.
locationInMesh      (-8.33 0.0 0.0); 
```

CfMesh also allows slightly different approaches, for example defining the `nCells` instead of `maxCellSize` or giving the bounding box size directly in the dictionary.


The next blocks are dedicated to **refinement inside the mesh**. Their role is to locally densify the mesh in specific areas. There are two ways to tell cfMesh how much it should be refined: providing cellSize (the maximum edge size of cells inside the area) or level (how many times the area should be refined from its current level). 

The first part called `objectRefinements` allows to place arbitrarily geometrical objects with refined cells. There are two basic types of objects: `box` and `sphere` which must have their origin, dimensions and `cellSize` or `level` provided.

```text
objectRefinements
{
    // name of the object, next inside the dictionary its parameters are provided.
    airfoilBox
    {
        type            box;
        centre         (4.25 0.0 0.0);
        lengthX        9.5;
        lengthY        3.0;
        lengthZ        1; // for 2D this parameter is negligible
        cellSize       0.01; // may be replaced with level
    }
    leadingEdgeSphere
    {
        type            sphere;
        centre          (0.0 0.0 0.0);
        radius          0.1;
        cellSize        0.0005;
    }

}
```

The next part called `localRefinement` looks at the local geometries defined inside the domain and refines them in a similar manner (by selecting the `cellSize` or `level`). The only difference is the presence of the `refinementThickness` key, which defines how far from the selected surface refining should be performed. In this situation the only surface available is the airfoil and it is important to stay consistent with patch names so cfMesh recognizes them correctly.

```text
localRefinement
{
    airfoil // name must be consistent with the previous patch definitions (eg. inside the FMS file)
    {
        cellSize 0.0015; // may be replaced with level
        refinementThickness 0.015;
    }
}
```

The last element is **boundary layers control**. Layers are very thin cells placed on some surface to catch the fluid behavior near the wall. They are not mandatory in every simulation to get reliable results (for example using the k-epsilon turbulence), but usually it is useful to have them, as they increase the accuracy without the need of increased global refinement. 

```text
boundaryLayers
{
    patchBoundaryLayers
    {
        airfoil
        {
            nLayers              18;
            thicknessRatio       1.2;
            maxFirstLayerThickness 5e-05;
            allowedDiscontinuity 1;
        }
    }
}
```

Having prepared the FMS file and a `meshDict` like in the example, you may try to launch the process by executing:

```powershell
cartesian2DMesh
```

#### FMS files

It is worth introducing some information about the FMS file format. Its main advantages are that it stores patches, subsets, and feature edges, making it a comprehensive source of geometric data for the meshing process. However, it is not widely used and generating it may be problematic. As mentioned before, for the airfoil cases automatic conversion is implemented, but here are a few steps for manual creation.

 - Creating the domain files in form of STL. In this scenario these are the airfoil, inlet, outlet, top and bottom surfaces. It is important to create them in ASCII format.
 - Combine them into a single domain STL file with separate solids. This may be achieved simply by copying them together into one file

```text
solid b'airfoil'
facet normal 0.000000 0.000000 -0.000001
  outer loop
    vertex 0.000000 0.000000 0.001000
    vertex 0.000625 0.003866 0.001000
    vertex 0.001050 0.005441 0.001000
  endloop
endfacet
...
solid b'inlet'
...
```

TODO

## Common meshing pathologies

### cfMesh z-span error

A common problem when trying to launch a semi-2D simulation is assuming the z-thickness of the mesh. In theory, it should have no effect, as there is only one cell in this direction. In practice, the thickness should scale with the minimum cell size of the mesh. Keeping the z-span about 10–100× the minimum cell size should be enough, but there is no universal rule. The problem may result in skewed cells created near the boundaries. Usually this will cause the simulation to crash in the first iteration, but sometimes it will be able to finish with visible numerical artifacts. 

Detailed information will show inside the `cartesian2DMesh` and `checkMesh` logs.

### STL resolution

Another problem that may result in skewed faces is the airfoil resolution. Crucial here is the leading edge which features high curvature. On one hand a well-resolved geometry here is required to obtain accurate results, on the other hand, too many samples may break the layer logic of cfMesh. A sign of this problem may be excessive meshing time.

### Boundary layers on leading edge

Even with a high-resolution STL, adding boundary layers to a leading edge with very high curvature can cause cells to "collapse" or become highly skewed. This happens when the requested layer thickness is too large relative to the local radius of curvature. There are two recommended solutions, one to add `objectRefinements` specifically around the leading edge, another to increase or add an additional instance in `localRefinement` at the airfoil patch.

## Boundary conditions

For a CFD simulation, every boundary must be precisely described. This means every relevant flow parameter needs an initial value and instructions on how to behave at the inlet, outlet, walls, etc. In OpenFOAM, this takes the form of text dictionary files dedicated to each field, placed inside the `0/` directory (time step 0).
Inside, all patches must have a proper boundary condition (for cfMesh they will be the same as defined in the `FMS` file). This means that the structure inside the files will look similar for almost all parameters. Below an example for `U`:

```C++
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
// specifies the dimension. In OpenFOAM: mass, length, time, temperature, quantity, current, luminous intensity
dimensions      [0 1 -1 0 0 0 0];

// only the initial value for the internal field
internalField   uniform (52.94117647058824 0 0);

// boundary condition types on all patches, names must correspond to patch definitions in mesh
boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform (52.94117647058824 0 0);
    }
    outlet
    {
        type            zeroGradient;
    }
    lowerWall
    {
        type            zeroGradient;
    }
    upperWall
    {
        type            zeroGradient;
    }
    front
    {
        type            symmetryPlane;
    }
    back
    {
        type            symmetryPlane;
    }
    airfoil
    {
        type            noSlip;
    }
}
```

Here are three important sections. First we have the **dimensions**. OpenFOAM uses 7 values with corresponding units: mass (kilograms), length (meters), time (seconds), temperature (kelvins), quantity (moles), current (amperes), luminous intensity (candelas). The list breaks down the above example.

```text
dimensions      [0 1 -1 0 0 0 0]; // kg^0, m^1, s^-1, K^0, mol^0, A^0, cd^0 -> the velocity unit m/s 
```

Similarly for pressure:
\[
[0\ 2\ -2\ 0\ 0\ 0\ 0]\ \Rightarrow\ \text{kg}^0\cdot\text{m}^2\cdot\text{s}^{-2}\cdot\text{K}^0\cdot\text{mol}^0\cdot\text{A}^0\cdot\text{cd}^0
\]
The resulting unit is \(\text{m}^2/\text{s}^2\), i.e. kinematic pressure \((p/\rho)\). In OpenFOAM, pressure is kinematic by default, do not confuse it with static pressure in Pa.

Next is the definition of the **internal field**. Every value must be initialized inside the simulation domain. There are several ways to do that, but in general initializations are uniform and it is recommended to try to make estimates according to the simulation assumptions. In many situations the field will converge to the accurate results from different internal fields, but having a good initialization may speed up the initial iterations, while having a bad one may lead to errors in the final result. Again let's take a look at the pressure field initialization. 

```text
internalField   uniform 0.0;        // for pressure (scalar)
internalField   uniform (50 0 0);   // for velocity (vector)
```

There are two important things to notice. First, as the velocity is a vector it is initialized in 3D space. Second, pressure is initialized as 0. In some situations, like in many airfoil simulations, the absolute pressure does not matter, rather the difference between two points. In such a situation the 0 relative pressure is often a clear approach.

Lastly the **boundary conditions**, which are the most important element inside the file. These are mathematical instructions on how the parameters should behave on the boundaries. This then determines processes inside the domain. Each patch must have a boundary specified, they should not contradict each other and should make physical sense. 

**Guidelines for setting boundary conditions:**
 - The flow must have a driving force: either specify velocity (Dirichlet) at inlet/outlet, or impose a pressure difference.
 - Avoid over-constraining: specify velocity *or* pressure at a given boundary, not both.
 - At solid walls, velocity must be zero (no-slip condition).
 - Ensure physical consistency: inlet conditions should match your Reynolds number and flow assumptions.

### Common boundary condition types

There are a lot of different boundary conditions in OpenFOAM, but there are two main groups of them:
 - Dirichlet - You tell the solver exactly what the value is at the wall (e.g., $U = 50\,m/s$).
 - Neumann - You tell the solver that the value at the wall is the same as the value in the cell next to it (the change is zero).
 
Here are some of the most popular implementations:

- **fixedValue**: Dirichlet condition that imposes a specified value on the patch (uniform or nonuniform). Typical for `U` at inlets or `T` on walls.
- **zeroGradient**: Neumann condition with zero normal gradient, the value is extrapolated from the interior. Common for `U` at outlets and `p` at inlets.
- **noSlip**: Velocity-only shorthand for fixedValue (0 0 0) on solid walls, enforces zero velocity at the wall.
- **symmetryPlane**: Symmetry/mirroring plane with zero normal components and zero normal gradients for tangential components, used for 2D front/back planes instead of `empty`, which is in some applications more problematic.
- **nutUSpaldingWallFunction**: Wall function for kinematic turbulent viscosity (nut) using the Spalding law-of-the-wall, applied on solid walls with wall-function RANS models.
- **kqRWallFunction**: Wall function for turbulence kinetic energy k at walls, enforces appropriate near-wall behavior and zero flux consistent with roughness model.
- **omegaWallFunction**: Wall function for specific dissipation rate ω at walls, sets ω based on y+, ν, and wall distance for k-ω/SST models.

## Numerical settings

The `system/` directory contains files that define how the continuous partial differential equations are discretized and solved numerically. The two most critical files are `fvSchemes` and `fvSolution`.

### fvSchemes: Discretization schemes

The `fvSchemes` dictionary specifies the numerical schemes used for discretization of the various terms in the Navier–Stokes equations. These choices directly affect accuracy and stability. The general rule is that first order schemes (like Gauss Upwind) are more stable, but less accurate, while second order schemes (like Gauss Linear) are more sensitive and tend to oscillate, but may give better results. 

```c++
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    // each solved value needs a specified divergence scheme
    default         none;
    div(phi,U)      Gauss upwind;
    div(phi,k)      Gauss upwind;
    div(phi,epsilon) Gauss upwind;
    div(phi,omega)   Gauss upwind;
    div(phi,nuTilda) Gauss upwind;
    div(phi,gammaInt) Gauss upwind;
    div(phi,ReThetat) Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

wallDist
{
    method          meshWave;
}
```

**Time derivatives (`ddtSchemes`):** Since airfoil simulations target steady-state results, the `default` is set to `steadyState`, meaning the transient term \(\partial/\partial t\) is omitted.

**Gradient schemes (`gradSchemes`):** Defines how gradients (\(\nabla\)) are calculated. The `Gauss linear` setting uses central differencing on the mesh, which is accurate for smooth fields.

**Divergence schemes (`divSchemes`):** Handle the convective terms (\(\mathbf{u}\cdot\nabla\)). The `Gauss upwind` scheme is used for velocity and turbulence properties to ensure stability by biasing calculations toward the upstream direction. The stress term uses `Gauss linear` for better accuracy.

**Laplacian schemes (`laplacianSchemes`):** Used for diffusion and viscous terms (\(\nabla^2\)). The `Gauss linear corrected` setting accounts for non-orthogonal mesh effects, improving accuracy when mesh angles deviate from 90°.

**Surface normal gradients (`snGradSchemes`):** Set to `corrected` to maintain accuracy on non-orthogonal meshes, particularly important near the curved leading edge.

### fvSolution: Solver and algorithm controls

The `fvSolution` file controls linear equation solvers, convergence tolerances, and the pressure–velocity coupling algorithm. Below a file example:

```c++
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
    // matrix solver specifications
    p
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-06;
        relTol          0.1;
    }
    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
    "(k|omega|epsilon|nuTilda)"
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }
    ReThetat
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0.3;
    }
    gammaInt
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.2;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 2;
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.7;
        k               0.6;
        omega           0.6;
        epsilon         0.7;
        nuTilda         0.7;
        ReThetat       0.3;
        gammaInt       0.2;
    }
}

cache
{
    grad(U);
}
```

#### Matrix solvers

Each field is assigned a specific solver and tolerance to balance accuracy and convergence speed.

 - **Pressure (`p`):** Uses the **GAMG** (Geometric-Agglomerated Algebraic Multi-Grid) solver with a `DIC` (Diagonal Incomplete Cholesky) smoother. This is highly efficient for pressure systems. The tight `tolerance` of `1e-06` ensures accurate pressure fields.
 
 - **Velocity (`U`):** Uses a `smoothSolver` with the `symGaussSeidel` smoother. The `tolerance` of `1e-05` is sufficient for momentum equations, `relTol` of `0` means residuals must meet the absolute tolerance.
 
 - **Turbulence parameters (`k`, `omega`, `epsilon`, `nuTilda`):** Use `smoothSolver` with Gauss–Seidel iterations to handle the tightly coupled nature of turbulence equations. The `relTol` of `0.1` allows faster convergence once residuals drop significantly.
 
 - **Transition parameters (`ReThetat`, `gammaInt`):** Use tight tolerances and lower relative tolerances to capture the precise moment and location of laminar–turbulent transition.

#### SIMPLE algorithm

The `SIMPLE` (Semi-Implicit Method for Pressure-Linked Equations) dictionary governs pressure–velocity coupling in steady-state cases.

 - **nNonOrthogonalCorrectors:** Set to `2` to improve convergence on non-orthogonal meshes. Each corrector iteration refines the pressure field to better account for mesh distortions near the airfoil's leading edge.

#### Relaxation factors

Because the Navier–Stokes equations are non-linear, solving them too aggressively between iterations can cause divergence. Relaxation factors "slow down" updates to maintain stability.

 - **Pressure (`p`):** Relaxed at `0.3`—heavily damped because pressure drives the entire flow solution.
 
 - **Velocity (`U`):** Relaxed at `0.7`—allows reasonable momentum updates while preventing instability.
 
 - **Turbulence kinetic energy and dissipation (`k`, `epsilon`, `omega`):** Relaxed at `0.6–0.7`—these equations are coupled to velocity and must be updated cautiously.
 
 - **Transition (`ReThetat`, `gammaInt`):** Relaxed at `0.3` and `0.2` respectively—very low to prevent oscillations during the detection of the transition point, ensuring smooth laminar–turbulent flow capture.

The `cache` directive pre-computes gradients of velocity to avoid redundant calculations each iteration.

## Simulation controls 

The `controlDict` file, located in the `system/` directory, is the central configuration file that governs the execution of the OpenFOAM solver. It defines the time-stepping parameters, data output frequency, and additional "on-the-fly" calculations known as function objects. Below a full example:

```c++
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}
// selected solver
application     simpleFoam;
// time control settings
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1500;
deltaT          1;
// write control settings
writeControl    timeStep;
writeInterval   250;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{
    // function objects for aerodynamic force coefficients
    forceCoeffs
    {
        type            forceCoeffs;
        libs            ("libforces.so");
        patches         (airfoil);
        rho         rhoInf;
        rhoInf          1.225;
        log             yes;
        writeControl    timeStep;
        writeInterval   5;
        CofR            (0.14195774447807372 0.012419693341541287 0);
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        pitchAxis       (0 0 1);
        magUInf         52.63157894736843;
        lRef            0.57;
        Aref            0.00011399999999999999;
        writeFields     no;
    }
}
```

### Time and Run Control

In OpenFOAM, even steady-state simulations (like those using `simpleFoam`) utilize a "time" variable. However, in steady-state, "time" actually represents the **iteration count**.

```text
application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1500;
deltaT          1;
```

* **application:** Specifies the solver to be used (e.g., `simpleFoam`).
* **startFrom / startTime:** Tells the solver whether to start from the beginning (`startTime`) or from the latest saved result (`latestTime`). This is useful for resuming a simulation that was stopped.
* **stopAt / endTime:** Defines the duration of the run. For airfoil simulations, you typically run for enough iterations (e.g., 1000–5000) for the residuals to drop and the coefficients ($C_l, C_d$) to stabilize.
* **deltaT:** For steady-state solvers, this is usually set to `1`, meaning the solver advances one iteration at a time.

### Data Management (Output)

CFD simulations can generate massive amounts of data. The `controlDict` allows you to manage how much of that data is kept.

```text
writeControl    timeStep;
writeInterval   250;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
```

* **writeControl / writeInterval:** Defines how often the solver saves the results (velocity, pressure, etc.) to the disk. For example, if `writeInterval` is 500, the solver saves every 500th iteration.
* **purgeWrite:** A vital setting to save disk space. It tells OpenFOAM to keep only a specific number of recent result folders, deleting older ones automatically.
* **writeFormat / writeCompression:** Using `binary` format and setting `writeCompression` to `on` can significantly reduce the file size of your results without losing numerical accuracy.

### Function Objects

Function objects are powerful utilities that allow you to perform calculations *during* the simulation rather than waiting for it to finish. In this repository, they are used primarily for aerodynamic monitoring.

#### Force Coefficients
To calculate lift and drag, a `forceCoeffs` function object is included. It integrates the pressure and skin friction over the airfoil surface patches and outputs the results to a text file in the `postProcessing/` directory.

**Key parameters required:**
- **patches:** The name of the airfoil surface (e.g., `airfoil`).
- **rhoInf:** The reference density. Must be consistent with simulation properties.
- **magUInf:** The free-stream velocity.  Must be consistent with simulation properties.
- **lRef / ARef:** The chord length and reference area (chord × span) used to normalize the forces into coefficients.
- **liftDir / dragDir:** The vectors defining the lift and drag directions relative to the coordinate system.

```c++
functions
{
    forces
    {
        type            forceCoeffs;
        libs            ("libforces.so");
        patches         (airfoil);
        rho             rhoInf;
        rhoInf          1.225;
        log             yes;
        writeControl    timeStep;
        writeInterval   5;
        CofR            (0.14195774447807372 0.012419693341541287 0);
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        pitchAxis       (0 0 1);
        magUInf         52.63157894736843;
        lRef            0.57;
        Aref            0.00011399999999999999;
        writeFields     no;
    }
}
```

#### Residuals Monitoring

Function objects are also used to track the convergence of the simulation in real-time. By monitoring the `postProcessing/residuals` file, you can observe how the error for each field ($p$, $\mathbf{u}$, $k$, $\omega$, etc.) evolves. 

* **Interpretation:** A declining residual indicates that the solver's guess for the flow field is becoming more accurate with each iteration.
* **Convergence Criteria:** If the residuals stop changing and reach a sufficiently low level—usually below $10^{-4}$ or $10^{-5}$—the simulation is considered successfully converged. At this point, the aerodynamic forces ($C_l, C_d$) should also show stable, constant values.
* **Troubleshooting:** If residuals are oscillating or increasing, it typically points to numerical instability, often caused by a poor-quality mesh near the leading edge or relaxation factors that are too high for the current flow conditions.

# OpenFOAM airfoil simulation

This chapter will be devoted to simulating airfoil using the repository. Short introduction, explanation about the contents.

## Repository structure

Diagram with short descriptions

## Workflow structure

Step by step typical simulation workflow description. Starting from preparation, meshing, simulating and postprocess. Mentioning the bash scripts that handle the parts on Ubuntu. Short mention about the outputs.

## Preparing case

Starting with the settings json that allow for global management. Introducing the python scripts that prepare the whole structure and files.

## Meshing

Mentioning the most important rules behind the meshing: bounding box size, where refinement for airfoils should be applied, short tutorial how to estimate cell size for different types of simulation (mostly depending on the turbulence models). How to look at the output, how to check for potential errors.

## Simulation

Preferred solvers, running parallel. What to look for to recognize potential problems.

## Postprocessing

Explaining the executed OpenFOAM commands. Listing and explaining the created plots. How to look at the results. Verifying residuals, results stability (coeffs plot in time), \(y^+\) and ParaView visualization. 

## Launching multiple cases

Demonstrating the code dedicated for multiple launches on different input jsons.

# Examples

## NACA4415 study

Launching the low Reynolds simulations with different angles of attack. Comparing them to the xfoil reference. Creating polar plots.

## Turbulence comparison

Launching the same case with different turbulence models.

## Demonstrating input json composition

Two same cases, where one has a minimum input json where all missing parameters are automatically filled with default values in the code and another where all this values are directly mentioned in the json.

# Summary

## Limitations

What is not supported here.

## Future development

Mentioning the features that are missing for now and potential areas of further development

## Useful

List some useful sites, material etc.
