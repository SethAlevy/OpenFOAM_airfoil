# OpenFOAM airfoil simulator

This repository is dedicated to perform 2D airfoil simulations using OpenFOAM (currently OpenFOAM v2406) and post-process. There are multiple functionalities to handle airfoils (generating, loading, exporting), create meshes, prepare boundary conditions and the whole case structure. Predefined environments and scripts allow to launch full analyses for several variants with post-processing and comparisons. 

# Getting started

To avoid setup difficulties and ensure repeatability in the performed analyses this repository provides predefined environments with only few prerequisites. Python environment is defined using Poetry, while for OpenFOAM a Docker image was created. 

The only truly requirement is Docker. It can be downloaded and installed using the [official page](https://www.docker.com/get-started/). Once the installation and configuration process is finished the image may be build. 

While inside the repository launch Docker Desktop and execute the following command in the terminal:

```powershell
docker build --progress=plain -t openfoam-airfoil .
```

This builds the image specified inside the [`Dockerfile`](Dockerfile). Starting from Ubuntu 22.04 all system dependencies, Python3.11 with Poetry and OpenFOAM v2406 with cf-mesh are installed on a container. The _--progres=plain_ parameter is optional and provides additional logs to track eventually build errors. All scripts inside the repository may be executed and used inside it. Usually it is useful to launch it with an attached local directory where the OpenFOAM cases and json files with simulation settings will be located. This allows to aces the files from the file explorer directly. This can be done with the following command:

```powershell
docker run -it -v "/input/dir/path:/app/case_dir" openfoam-airfoil
```

Now the repository and the case directory may be found under the _/app_ path. For example the main cases running script will be under the path: /app/src/simulation/run/run_case.sh. One may launch it using one of the example inputs by the following commands:

```powershell
cd /app/src/simulation/run/
bash ./run_case.sh --working-dir /app/case_dir --setup-file /app/examples/naca4415_study/aoa_5.json --case-name aoa_5
```

This example launches a full workflow for a generated NACA4415 airfoil with an angle of attack of 5 deg and Reynolds number of 1000000 set tu the boundary condition.  

The scripts are prepared to utilize Python scripts by the Poetry environment and it can be used also manually inside the container, but sometimes it is useful to set it up locally. A concise installation guide is available in its [documentation](https://python-poetry.org/docs/). One of the key project files for Poetry is [`pyproject.toml`](pyproject.toml) in the repository root. It defines required libraries and their versions. To ensure compatibility and proper operation, versions are pinned.

To build the environment, execute:

```powershell
poetry install
```

Now it may be used through the code editor for notebooks and scripts or through the terminal:

```powershell
poetry run python src/simulation/preparation/prepare_case.py --working-dir . --setup-file examples/naca4415_study/aoa_5.json --case-name case_name
```


# Airfoil

In general an airfoil is an aerodynamic shape specially designed to generate lift (due to pressure difference on its lower and upper surface) and minimum drag while moving through the air. They can be found wings or propellers as their cross sections.

## Airfoil terminology

There is a set of common parameters that are used to describe and understand the geometry and aerodynamic behavior of most airfoils. The most important terms are listed below:

 - Leading edge - the foremost point of the airfoil, that encounters the incoming flow.
 - Trailing edge - the rearmost point of the airfoil, where the upper and lower surfaces meet and the flow leaves the airfoil.
 - Chord - a straight line joining the leading edge and trailing edge of the airfoil. The chord length is  one of the most important geometric parameters, commonly used as a reference dimension in aerodynamic equations.
 - Camber - a curve equivalent to the the geometrical centerline. LLocated midway between the upper and lower surfaces describing the airfoils curvature.
 - Thickness - the distance between the upper and lower surfaces. Thickness is most commonly measured perpendicular to the camber line, though it may also be defined perpendicular to the chord line. The maximum thickness is a frequently used geometric parameter.
 - Upper and lower surfaces - create the airfoil outer geometry connecting the leading and trailing edges. Typically the upper surface encounters lower static pressure than the lower and sometimes they are mentioned as suction and pressure surfaces. The pressure difference between them generates lift. 
 - Aerodynamic center - a point along the chord where the the pitching moment stays unchanged regardless the change of angle of attack and fluid speed.
 - Center of pressure - a point where the average pressure force is considered to act. Its location varies with the changing conditions.
 - Angle of attack - the angle between the relative wind vector and the chord line. 
 - Lift - the component of the aerodynamic force acting perpendicular to the incoming flow direction.
 - Drag - the component of the aerodynamic force acting parallel to the incoming flow direction.
 
Place for plot

## NACA airfoils 

The National Advisory Committee for Aeronautics (NACA) developed and tested a series of airfoils in the first half of the 20th century. They were designed to provide systematic, well-defined geometries that could be easily reproduced, analyzed, and compared in experiments and calculations. Their most important feature is the description of key geometrical parameters through mathematical equations. Main parameters may be determined through the designations. Although there are modern, more capable airfoils developed the NACA-series are widely used for education, aerodynamic studies, numerical simulations, validations and experiments. There are two main groups of NACA airfoils.

## NACA 4-digits - the simplest series where the geometry is described by three parameters through four digits in form of:

**NACA MPXX**

where: 
 - **M** is the maximum camber as a percentage of the chord, 
 - **P** is the position of maximum camber as tenths of the chord,
 - **XX** is the maximum thickness as a percentage of the chord. 
 
For example NACA 4415 (4-4-15) has a maximum camber of 4% of the chord, located at 40% of the chord line from the leading edge and maximum thickness of 15% of the chord. Another good example is the NACA 0012 (0-0-12) which has a zero camber (which means the airfoil is symmetric and does not produce lift at 0 deg angle of attack) and maximum thickness of 12% of the chord. The following equations allow to calculate the airfoils geometry: 

### Normalized parameters

Using a normalized chordwise coordinate \( x \in [0,1] \):

\[
m = \frac{M}{100}, \quad
p = \frac{P}{10}, \quad
t = \frac{XX}{100}
\]

---

### Mean camber line

\[
y_c(x) =
\begin{cases}
\frac{m}{p^2}(2px - x^2), & x \le p \\
\frac{m}{(1-p)^2}\left[(1 - 2p) + 2px - x^2\right], & x > p
\end{cases}
\]

---

### Thickness distribution

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

### Surface inclination

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

### Upper and lower surface coordinates

\[
x_u = x - y_t \sin\theta, \quad
y_u = y_c + y_t \cos\theta
\]

\[
x_l = x + y_t \sin\theta, \quad
y_l = y_c - y_t \cos\theta
\]

---

## NACA 5-digits - slightly more complex in geometry and description than the 4-digits. Their naming convention is in form of:

**NACA LPQXX** 

where: 
 - **L** after multiplying by 3/20 is the design ideal lift coefficient, 
 - **P** after multiplying by 0.05 is the location of the maximum camber as a tenth of the chord, 
 - **Q** identifies the type of the camber (0 if simple, 1 if reflexed), 
 - **XX** is the maximum thickness as a percentage of the chord. 
 
For example NACA 23012 (2-3-0-12) has a design lift coefficient of 0.3, maximum camber located at 15% of the chord, maximum thickness of 12% of the chord and a simple camber.

### Normalized parameters

\[
c_l = 0.15L, \quad
p = 0.05P, \quad
t = \frac{XX}{100}
\]

The thickness distribution is identical to the NACA 4-digit series.

---

### Mean camber line (normal camber, Q = 0)

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

---

### Mean camber line (reflex camber, Q = 1)

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

---

### Surface coordinates

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

---


## Other NACA airfoil series

Additional NACA airfoil families, such as the 6-digit and 16-digit series, were developed for more specialized aerodynamic performance, particularly with respect to pressure distribution and laminar flow control. These series use more complex definitions and are not supported by the current implementation, and therefore are not discussed further.

## UIUC database

Maintained by the University of Illinois the site is a widely used source of information about airfoils and other related to aerodynamic resources. At this moment it provides more than 1600 geometries that may be downloaded in form of a dat file. This form allows automatic implementation into the repository knowing only the exact designation of the desired geometry. 

## Repository airfoil tools

The repository has as set of implemented tools that aim to make the work with airfoils easier. There are two options to get a geometry: 

 - generating a NACA 4-digit or NACA 5-digit airfoil
 - import an airfoil from UIUC database by designation. 

