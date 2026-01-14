def _patch(name, setup, default):
    return name if name is not None else setup.get(default, default)


def U_bc(
    velocity_str,
    inlet_patch=None,
    outlet_patch=None,
    lower_wall_patch=None,
    upper_wall_patch=None,
    front_patch=None,
    back_patch=None,
    airfoil_patch=None,
    setup=None
):
    setup = setup or {}
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform ({velocity_str});

boundaryField
{{
    {_patch(inlet_patch, setup, "inlet")}
    {{
        type            fixedValue;
        value           uniform ({velocity_str});
    }}
    {_patch(outlet_patch, setup, "outlet")}
    {{
        type            zeroGradient;
    }}
    {_patch(lower_wall_patch, setup, "lowerWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(upper_wall_patch, setup, "upperWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(front_patch, setup, "front")}
    {{
        type            symmetryPlane;
    }}
    {_patch(back_patch, setup, "back")}
    {{
        type            symmetryPlane;
    }}
    {_patch(airfoil_patch, setup, "airfoil")}
    {{
        type            noSlip;
    }}
}}
"""


def p_bc(
    pressure_value,
    inlet_patch=None,
    outlet_patch=None,
    lower_wall_patch=None,
    upper_wall_patch=None,
    front_patch=None,
    back_patch=None,
    airfoil_patch=None,
    setup=None
):
    setup = setup or {}
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform {pressure_value};

boundaryField
{{
    {_patch(inlet_patch, setup, "inlet")}
    {{
        type            zeroGradient;
    }}
    {_patch(outlet_patch, setup, "outlet")}
    {{
        type            fixedValue;
        value           uniform {pressure_value};
    }}
    {_patch(lower_wall_patch, setup, "lowerWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(upper_wall_patch, setup, "upperWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(front_patch, setup, "front")}
    {{
        type            symmetryPlane;
    }}
    {_patch(back_patch, setup, "back")}
    {{
        type            symmetryPlane;
    }}
    {_patch(airfoil_patch, setup, "airfoil")}
    {{
        type            zeroGradient;
    }}
}}
"""


def k_bc(
    k,
    inlet_patch=None,
    outlet_patch=None,
    lower_wall_patch=None,
    upper_wall_patch=None,
    front_patch=None,
    back_patch=None,
    airfoil_patch=None,
    setup=None
):
    setup = setup or {}
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform {k};

boundaryField
{{
    {_patch(inlet_patch, setup, "inlet")}
    {{
        type            fixedValue;
        value           uniform {k};
    }}
    {_patch(outlet_patch, setup, "outlet")}
    {{
        type            zeroGradient;
    }}
    {_patch(lower_wall_patch, setup, "lowerWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(upper_wall_patch, setup, "upperWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(front_patch, setup, "front")}
    {{
        type            symmetryPlane;
    }}
    {_patch(back_patch, setup, "back")}
    {{
        type            symmetryPlane;
    }}
    {_patch(airfoil_patch, setup, "airfoil")}
    {{
        type            kqRWallFunction;
        value           uniform 0;
    }}
}}
"""


def omega_bc(
    omega,
    inlet_patch=None,
    outlet_patch=None,
    lower_wall_patch=None,
    upper_wall_patch=None,
    front_patch=None,
    back_patch=None,
    airfoil_patch=None,
    setup=None
):
    setup = setup or {}
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      omega;
}}

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform {omega};

boundaryField
{{
    {_patch(inlet_patch, setup, "inlet")}
    {{
        type            fixedValue;
        value           uniform {omega};
    }}
    {_patch(outlet_patch, setup, "outlet")}
    {{
        type            zeroGradient;
    }}
    {_patch(lower_wall_patch, setup, "lowerWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(upper_wall_patch, setup, "upperWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(front_patch, setup, "front")}
    {{
        type            symmetryPlane;
    }}
    {_patch(back_patch, setup, "back")}
    {{
        type            symmetryPlane;
    }}
    {_patch(airfoil_patch, setup, "airfoil")}
    {{
        type            omegaWallFunction;
        value           uniform 0;
    }}
}}
"""


def epsilon_bc(
    epsilon,
    inlet_patch=None,
    outlet_patch=None,
    lower_wall_patch=None,
    upper_wall_patch=None,
    front_patch=None,
    back_patch=None,
    airfoil_patch=None,
    setup=None
):
    setup = setup or {}
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      epsilon;
}}

dimensions      [0 2 -3 0 0 0 0];

internalField   uniform {epsilon};

boundaryField
{{
    {_patch(inlet_patch, setup, "inlet")}
    {{
        type            fixedValue;
        value           uniform {epsilon};
    }}
    {_patch(outlet_patch, setup, "outlet")}
    {{
        type            zeroGradient;
    }}
    {_patch(lower_wall_patch, setup, "lowerWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(upper_wall_patch, setup, "upperWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(front_patch, setup, "front")}
    {{
        type            symmetryPlane;
        value           uniform {epsilon};
    }}
    {_patch(back_patch, setup, "back")}
    {{
        type            symmetryPlane;
        value           uniform {epsilon};
    }}
    {_patch(airfoil_patch, setup, "airfoil")}
    {{
        type            fixedValue;
        value           uniform 0;
    }}
}}
"""


def nut_bc(
    inlet_patch=None,
    outlet_patch=None,
    lower_wall_patch=None,
    upper_wall_patch=None,
    front_patch=None,
    back_patch=None,
    airfoil_patch=None,
    setup=None
):
    setup = setup or {}
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}}

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
    {_patch(inlet_patch, setup, "inlet")}
    {{
        type            zeroGradient;
    }}
    {_patch(outlet_patch, setup, "outlet")}
    {{
        type            zeroGradient;
    }}
    {_patch(lower_wall_patch, setup, "lowerWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(upper_wall_patch, setup, "upperWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(front_patch, setup, "front")}
    {{
        type            symmetryPlane;
    }}
    {_patch(back_patch, setup, "back")}
    {{
        type            symmetryPlane;
    }}
    {_patch(airfoil_patch, setup, "airfoil")}
    {{
        type            nutUSpaldingWallFunction;
        value           uniform 0;
    }}
}}
"""


def gammaInt_bc(
    gammaInt=1e-5,
    inlet_patch=None,
    outlet_patch=None,
    lower_wall_patch=None,
    upper_wall_patch=None,
    front_patch=None,
    back_patch=None,
    airfoil_patch=None,
    setup=None
):
    setup = setup or {}
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      gammaInt;
}}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform {gammaInt};

boundaryField
{{
    {_patch(inlet_patch, setup, "inlet")}
    {{
        type            fixedValue;
        value           uniform {gammaInt};
    }}
    {_patch(outlet_patch, setup, "outlet")}
    {{
        type            zeroGradient;
    }}
    {_patch(lower_wall_patch, setup, "lowerWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(upper_wall_patch, setup, "upperWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(front_patch, setup, "front")}
    {{
        type            symmetryPlane;
    }}
    {_patch(back_patch, setup, "back")}
    {{
        type            symmetryPlane;
    }}
    {_patch(airfoil_patch, setup, "airfoil")}
    {{
        type            zeroGradient;
    }}
}}
"""


def retheta_bc(
    reteta=1000,
    inlet_patch=None,
    outlet_patch=None,
    lower_wall_patch=None,
    upper_wall_patch=None,
    front_patch=None,
    back_patch=None,
    airfoil_patch=None,
    setup=None
):
    setup = setup or {}
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      ReTheta;
}}

dimensions      [0 0 0 0 0 0 0];

internalField   uniform {reteta};

boundaryField
{{
    {_patch(inlet_patch, setup, "inlet")}
    {{
        type            fixedValue;
        value           uniform {reteta};
    }}
    {_patch(outlet_patch, setup, "outlet")}
    {{
        type            zeroGradient;
    }}
    {_patch(lower_wall_patch, setup, "lowerWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(upper_wall_patch, setup, "upperWall")}
    {{
        type            zeroGradient;
    }}
    {_patch(front_patch, setup, "front")}
    {{
        type            symmetryPlane;
    }}
    {_patch(back_patch, setup, "back")}
    {{
        type            symmetryPlane;
    }}
    {_patch(airfoil_patch, setup, "airfoil")}
    {{
        type            fixedValue;
        value           uniform 1000;
    }}
}}
"""
