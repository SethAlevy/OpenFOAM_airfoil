"""
SnappyHexMesh generation utilities for OpenFOAM.
"""
from pathlib import Path
from airfoil.airfoil import Airfoil
from templates.initial_settings_template import Settings
from simulation.openfoam.block_mesh import get_bounding_box


def snappy_hex_mesh_dict(
        airfoil: Airfoil,
        setup: Settings,
        output_path: Path
) -> None:
    """
    Generate snappyHexMeshDict for OpenFOAM.

    Args:
        airfoil (Airfoil): Airfoil object containing geometry.
        setup (Settings): Settings object containing mesh parameters.
        output_path (Path): Path to save the snappyHexMeshDict file.
    """
    mesh_settings = setup.mesh_settings
    snappy = mesh_settings.get("SnappyHexMesh", {})

    x_min, x_max, y_min, y_max, z_min, z_max = get_bounding_box(
        airfoil, mesh_settings.get("BoundingBox", {}))

    location_in_mesh = (x_min + 0.1 * airfoil.chord, 0, 0)

    refinement_min = snappy.get("RefinementMinLevel", 2)
    refinement_max = snappy.get("RefinementMaxLevel", 4)

    content = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}

castellatedMesh true;
snap            true;
addLayers       false;

geometry
{{
    airfoil.stl
    {{
        type triSurfaceMesh;
        name airfoil;
    }}
}};

castellatedMeshControls
{{
    maxLocalCells 100000;
    maxGlobalCells 2000000;
    minRefinementCells 0;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 3;
    
    features
    (
    );
    
    refinementSurfaces
    {{
        airfoil
        {{
            level ({refinement_min} {refinement_max});
        }}
    }}
    
    resolveFeatureAngle 30;
    
    refinementRegions
    {{
    }}
    
    locationInMesh ({location_in_mesh[0]} {location_in_mesh[1]} {location_in_mesh[2]});
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 30;
    nRelaxIter 5;
}}

addLayersControls
{{
    relativeSizes true;
    layers
    {{
    }}
    expansionRatio 1.0;
    finalLayerThickness 0.3;
    minThickness 0.1;
    nGrow 0;
    featureAngle 30;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}}

meshQualityControls
{{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-13;
    minTetQuality 1e-30;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.02;
    minVolRatio 0.01;
    minTriangleTwist -1;
    
    nSmoothScale 4;
    errorReduction 0.75;
}}

debug 0;
mergeTolerance 1e-6;
"""

    with open(output_path / "snappyHexMeshDict", "w") as f:
        f.write(content)

