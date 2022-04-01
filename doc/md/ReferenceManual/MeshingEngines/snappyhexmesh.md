@page snappyhexmesh_ref SnappyHexMesh

[TOC]

\section desc Description

<em>snappyHexMesh</em> is a meshing engine that has been developed specifically for quality CFD hex(-dominant) mesh 
generation. It is included with <em>OpenFOAM</em> and requires that the environment has been properly set. On a default
installation, this can be done by executing the following command in the working terminal:

    . /usr/lib/openfoam/openfoam2006/etc/bashrc

Mesh generation with <em>snappyHexMesh</em> involves three steps:
1. <strong>Castellating</strong>:  removes any cell beyond the meshing region of interest set by the user.
2. <strong>Snapping</strong>: moves the jagged edges of the cell along the boundary so that they are coincident with the boundary.
3. <strong>Layering</strong>: adds layers in the boundary region to ensure higher-quality cells in the transition region.

The tutorial [available here](https://cfd.direct/openfoam/user-guide/v9-snappyHexMesh/) provides an excellent visualization
of how these three steps are applied to a 2D shape.

<em>snappyHexMesh</em> requires a purely hexahedral base mesh, which can be created with another <em>Promesh</em> tool 
such as <em>blockMesh</em>.

\section ex Examples

Generic code for <em>snappyHexMesh</em>:
```json
{
  "Program Type": "Mesh Generation",
  "Mesh File Options": {
    "Input Geometry File": "geometry.stl"
    "Output Mesh File": "meshed_geometry.vtu"
  },
  "Mesh Generation Options": {
    "Mesh Generation Engine": "snappyHexMesh",
    "snappyHexMesh Parameters": {
      "Castellated Mesh": true,
      "Snapping": true, 
      "Layer Addition": false,
      "mergeTolerance": 1e-06,
      "Geometry Definition": {
        "Enable Multi Patches": true | false
        "Input Patch Name": "patch_name"
        "Custom Patches": [
          {
            "Custom Patch Name": "box_patch_name",
            "Searchable Shape": "searchableBox",
            "minimum bound": [-0.5, -0.5, -0.5],
            "maximmum bound": [0.5, 0.5, 0.5]
          },
          {
            "Custom Patch Name": "cylinder_patch_name",
            "Searchable Shape": "searchableCylinder",
            "Radius": 0.5,
            "Axis Point 1": [-0.5, 0.0, 0.0],
            "Axis Point 2": [0.5, 0.0, 0.0]
          },
          {
            "Custom Patch Name": "sphere_patch_name",
            "Searchable Shape": "searchableSphere",
            "Radius": 0.5,
            "Center": [0.0, 0.0, 0.0]
          }
        ]
      },
      "Castellated Mesh Controls": {
        "CellZones": false,
        "maxLocalCells": 2000000,
        "maxGlobalCells": 4000000,
        "minRefCells": 0,
        "General GapLevelIncrement": 1,
        "nCellsBetweenLevels": 3,
        "surfaceRefinementLvlMin": 0,
        "surfaceRefinementLvlMax": 0,
        "resolveFeatureAngle": 60,
        "gapLevelIncrement": 1,
        "planarAngle": 30,
        "locationInMesh": [0.0, 0.0, 0.0],
        "allowFreeStandingZoneFaces": true,
        "Feature File": {
          "File Name": "geometry_feature_edge_file.eMesh",
          "MinLevel": 1,
          "MaxLevel": 3
        },
        "SurfaceRefinementRegions": {
          "Patch Name": "patch_name_from_stl",
          "Patch Type": "patch" | "wall",
          "MinLevel": 1,
          "MaxLevel": 3
        },
        "GeomRefinementRegions": {
          "Patch Name": "custom_patch_name_from_geometry_def",
          "Mode": "inside" | "outside",
          "MinLevel": 1,
          "MaxLevel": 3
        }
      },
      "Snapping Controls": {
        "nSmoothPatch": 4,
        "tolerance": 0.5,
        "snapSolveIter": 200,
        "snapRelaxIter": 6,
        "nFeatureSnapIter": 10,
        "implicitFeatureSnap": false,
        "explicitFeatureSnap": false,
        "multiRegionFeatureSnap": false
      },
      "Mesh Layer Controls": {
        "relativeSizes": true,
        "expansionRatio": 1.3,
        "finalLayerThickness": 1,
        "minThickness": 0.1,
        "firstLayerThickness": 0.1,
        "thickness": 0.2,
        "nGrow": 0,
        "featureAngle": 30,
        "nRelaxIter": 3,
        "nSmoothSurfaceNormals": 1,
        "nSmoothNormals": 3,
        "nSmoothThickness": 2,
        "maxFaceThicknessRatio":  0.5,
        "maxThicknessToMedialRatio": 1,
        "minMedialAxisAngle": 90,
        "nBufferCellsNoExtrude": 0,
        "nLayerIter": 50,
        "nRelaxedIter": 20,
        "slipFeatureAngle": 30,
        "nMedialAxisIter": 1,
        "nSmoothDisplacement": 1,
        "Layers": [
          {
            "Patch Name": "custom_patch_combo_1",
            "nSurfaceLayers": 1,
            "minThickness": 0.1,
            "expansionRatio": 1,
            "finalLayerThickness": 1
          },
          {
            "Patch Name": "custom_patch_combo_2",
            "nSurfaceLayers": 1,
            "minThickness": 0.1,
            "expansionRatio": 1,
            "firstLayerThickness": 0.1
          },
          {
            "Patch Name": "custom_patch_combo_3",
            "nSurfaceLayers": 1,
            "minThickness": 0.1,
            "thickness": 1,
            "finalLayerThickness": 1
          },
          {
            "Patch Name": "custom_patch_combo_4",
            "nSurfaceLayers": 1,
            "minThickness": 0.1,
            "thickness": 1,
            "firstLayerThickness": 0.1
          },
          {
            "Patch Name": "custom_patch_combo_5",
            "nSurfaceLayers": 1,
            "minThickness": 0.1,
            "expansionRatio": 1,
            "thickness": 1
          }
        ]        
      },
      "Mesh Quality Controls": {
        "maxNonOrtho": 65,
        "maxBoundarySkewness": 20,
        "maxInternalSkewness": 4,
        "maxConcave": 80,
        "minVol": 1e-13,
        "minTetQuality": 1e-15,
        "minArea": -1,
        "minTwist": 0.02,
        "minFaceWeight": 0.05,
        "minVolRatio": 0.01,
        "minDeterminant": 0.001,
        "minTriangleTwist": -1,
        "qcnSmoothScale": 5,
        "errorReduction": 0.75
      }
    }
  }
}
```

\section params Parameters

\subsection general_params snappyHexMesh Parameters

- <strong>`Castellated Mesh`</strong>:  If `true`, enables castellated mesh option.  <strong>Required</strong> 
- <strong>`Snapping`</strong>:  If `true`, enables surface snapping option. <strong>Required</strong> 
- <strong>`Layer Addition`</strong>:  If `true`, enables adding layers. <strong>Required</strong>
- <strong>`mergeTolerance`</strong>:  If specified, specifies a merge tolerance as a fraction of the bounding box of the initial mesh.  Default is 1e-06 

\subsection geom_params Geometry Definition

- <strong>`Enable Multi Patches`</strong>:  Should be enabled if the geometry file has more than one patch.  <strong>Required</strong> 
- <strong>`Input Patch Name`</strong>:  If specified, will be assigned to the entire surface of the provided geometry.  Optional; should only be specified if entire .stl is one solid
- <strong>`Custom Patches`</strong>:  Array of custom patches defined with primitive shapes.  Optional 
    - <strong>`Custom Patch Name`</strong>:  User provided name given to the custom patch  
    - <strong>`Searchable Shape`</strong>:  One of `"searchable<SHAPE>"`, where `<SHAPE>` is `Box`, `Cylinder`, or `Sphere`.  
    - <strong>`minimum bound`</strong>:  An array of coordinates specifying the minimum corner of a `searchableBox`.  For example, [-35, -35, -35] 
    - <strong>`maximum bound`</strong>:  An array of coordinates specifying the maximum corner of a `searchableBox` (the opposite diagonal from the `minimum bound`).  For example, [-50, 35, 52] 
    - <strong>`Axis Point 1`</strong>:  An array of coordinates specifying the center of the base of a `searchableCylinder`.  For example, [0, 0, 0] 
    - <strong>`Axis Point 2`</strong>:  An array of coordinates specifying the center of the top of a `searchableCylinder`. For example, [0, 0, 35] 
    - <strong>`Radius`</strong>:  The radius of a `searchableCylinder` or a `searchableSphere`.  
    - <strong>`Center`</strong>:  An array of coordinates specifying the center of a `searchableSphere`.  

\subsection cast_mesh Castellated Mesh Controls

Refinement in <em>snappyHexMesh</em> happens via splitting cells in half. For a 2 cm background mesh, the first level of 
refinement will refine it to a 1 cm mesh. The second level of refinement will further refine it to a 0.5 cm mesh. 

Each level of refinement splits cells in half in each dimension. The general formula is \f$N(2^d)^R\f$, where \f$N\f$ is
the number of cells being split, \f$d\f$ is the number of dimensions (2 or 3), and \f$R\f$ is the levels of refinement
being applied. For example, in a 2D cell, the first level of refinement splits it into \f$(2^2)^1 = 4\f$ cells. The 
fifth level of refinement would produce \f$(2^2)^5 = 1024 \f$ cells, each \f$2^{(-5)}\f$ the size. Similarly, for a 3D 
cell, the first level of refinement splits it into \f$(2^3)^1 = 8\f$ cells. The fifth level of refinement would produce 
\f$(2^3)^5 = 32768\f$ cells.

- <strong>`CellZones`</strong>:  If `true`, indicates that there will be multiple disconnected regions within the mesh.  <strong>Required</strong> 
- <strong>`maxLocalCells`</strong>:  Specifies the maximum number of cells per processor.  Optional; default is 2000000 
- <strong>`maxGlobalCells`</strong>:  Specifies the global maximum number of cells.  Optional; default is 4000000 
- <strong>`minRefCells`</strong>:  Specifies the minimum number of refinement cells.  Optional; default is 0 
- <strong>`GeneralGapLevelIncrement`</strong>:    Optional; default is 1 
- <strong>`nCellsBetweenLevels`</strong>:   Sets the number of buffer layers between different levels of refinement.  Optional; default is 3 
- <strong>`surfaceRefinementLvlMin`</strong>:  Sets the minimum level of surface refinement.  Optional; default is 0 
- <strong>`surfaceRefinementLvlMax`</strong>:  Sets the maximum level of surface refinement.  Optional; default is 0 
- <strong>`resolveFeatureAngle`</strong>:  Sets the maximum intersection angle a cell can see---any higher, and the maximum level of refinement will be applied. Note that this is the inverse of the corresponding feature in <em>cfMesh</em>.  Optional; default is 60 degrees 
- <strong>`gapLevelIncrement`</strong>:    Optional; default is 1 
- <strong>`planarAngle`</strong>:    Optional; default is 30 degrees 
- <strong>`locationInMesh`</strong>:  Specifies where meshing will occur. If this point falls within a volume, that inner volume will be meshed. If it falls outside the volume, the outer volume will be meshed.  Optional; default is [0, 0, 0] 
- <strong>`allowFreeStandingZoneFaces`</strong>:  Should be set to `true` if the `faceZones` (specified in the `SurfaceRefinementRegions` section) can be free-standing rather than restricted to the boundary of the corresponding `cellZones`.  Optional; default is `true` 
- <strong>`Feature File`</strong>:    Optional 
    - <strong>`File Name`</strong>:  Name of the geometry feature edge file (.eMesh format), if it exists.  <strong>Required</strong> 
    - <strong>`MinLevel`</strong>:  Sets the minimum level of refinement at feature edges.  Optional; default is 1
    - <strong>`MaxLevel`</strong>:  Sets the maximum level of refinement at feature edges.  <strong>Required</strong> 
- <strong>`SurfaceRefinementRegions`</strong>:    Optional 
    - <strong>`Patch Name`</strong>:  Name (taken from the stl file) of the patch.  <strong>Required</strong> 
    - <strong>`Patch Type`</strong>:  Type of patch; can be either `"patch"` or `"wall"`.  Optional; default is `"NO"` 
    - <strong>`MinLevel`</strong>:  Sets the minimum level of refinement for the surface patch.  Optional; default is 1 
    - <strong>`MaxLevel`</strong>:  Sets the maximum level of refinement for the surface patch.  <strong>Required</strong> 
- <strong>`GeomRefinementRegions`</strong>:    Optional 
    - <strong>`Patch Name`</strong>:  One of the custom patch names defined in the `Geometry Definition` section.  
    - <strong>`Mode`</strong>:  Indicates whether refinement should take place `inside` or `outside` the patch.  
    - <strong>`MinLevel`</strong>:  Sets the minimum level of refinement for the geometry.  
    - <strong>`MaxLevel`</strong>:  Sets the maximum level of refinement for the geometry.  
    
\subsection snapping Snapping Controls

- <strong>`nSmoothPatch`</strong>:  Specifies the number of smoothing iterations to conform the mesh to the geometry surface.  Optional; default is 4 
- <strong>`tolerance`</strong>:  Sets how far the algorithm will search for a point to snap to. Specifically, the ratio of the distance between a point and the surface feature to the local maximum edge length.  Optional; default is 0.5 
- <strong>`snapSolveIter`</strong>:  Number of iterations for the snapping algorithm solver.  Optional; default is 200 
- <strong>`snapRelaxIter`</strong>:  Number of iterations for the snapping algorithm relaxation.   Optional; default is 6 
- <strong>`nFeatureSnapIter`</strong>:  Number of iterations for feature snapping  Optional; default is 10 
- <strong>`implicitFeatureSnap`</strong>:  If `true`, enables snapping based on sharp edges in the geometry.  Optional; default is `false` 
- <strong>`explicitFeatureSnap`</strong>:  If `true`, enables snapping based on feature edges defined in the `Feature File`, specified in the `Castellated Mesh Controls` section.  Optional; default is `false` 
- <strong>`multiRegionFeatureSnap`</strong>:  If `true` (and `explicitFeatureSnap` is enabled), features between multiple surfaces will be captured.  Optional; default is `false` 

\subsection mesh_layer Mesh Layer Controls

- <strong>`relativeSizes`</strong>:  Enables relative sizes option during layer addition. Allows the thickness of layers to be specified relative to the undistorted cell size rather than in absolute terms.  Optional; default is `true` 
- <strong>`expansionRatio`</strong>:  Expansion ratio for layer addition.  Optional; default is 1.3 
- <strong>`finalLayerThickness`</strong>:  Thickness (relative or absolute) of the layer furthest from the wall.  Optional; default is 1. 
- <strong>`minThickness`</strong>:  Minimum overall thickness of the layers.  Optional; default is 0.1 
- <strong>`firstLayerThickness`</strong>:  Thickness (relative or absolute) of the layer closest to the wall.  Optional; default is -1. 
- <strong>`thickness`</strong>:  Overall thickness of the cells. Ignored if `finalLayerThickness` or `firstLayerThickness` are specified.  Optional; default is -1. 
- <strong>`nGrow`</strong>:  Sets the number of layers that should not be grown if points are not extruded. Used to delay layer growth close to features.  Optional; default is 0 
- <strong>`featureAngle`</strong>:  The angle in degrees above which surfaces are not extruded for layer addition.  Optional; default is 30 degrees. If layers at corner edges are desired, set to 180 
- <strong>`nRelaxIter`</strong>:  Number of relaxation steps.  Optional; default is 3 
- <strong>`nSmoothSurfaceNormals`</strong>:  Number of patch normal smoothing operations.  Optional; default is 1 
- <strong>`nSmoothNormals`</strong>:  Number of smoothing iterations of interior mesh movement directions.  Optional; default is 3 
- <strong>`nSmoothThickness`</strong>:  Smooth layer thickness over the surface patches.  Optional; default is 2 
- <strong>`maxFaceThicknessRatio`</strong>:  Sets a maximum face to thickness ratio to stop layer growth on highly warped cells.  Optional; default is 0.5 
- <strong>`maxThicknessToMedialRatio`</strong>:  Sets a maximum thickness to medial distance ratio to reduce layer growth, typically in narrow cavities.  Optional; default is 1 
- <strong>`minMedialAxisAngle`</strong>:  Minimum angle to select the medial axis points.  Optional; default is 90 degrees 
- <strong>`nBufferCellsNoExtrude`</strong>:  Sets the number of cells in a buffer region to gradually step down the number of layers.  Optional; default is 0 
- <strong>`nLayerIter`</strong>:  Maximum number of layer addition iterations.  Optional; default is 50 
- <strong>`nRelaxedIter`</strong>:  Number of iterations after which relaxed mesh quality controls are used.  Optional; default is 20 
- <strong>`slipFeatureAngle`</strong>:  Allows the sliding of points on the patch without the layer grown if angle to the patch extrusion direction is larger than the one specified.  Optional; default is 30 degrees 
- <strong>`nMedialAxisIter`</strong>:  Limits number of steps walking away from the surface.  Optional; default is -1 
- <strong>`nSmoothDisplacement`</strong>:  Smoothing displacement of the mesh after the medial axis is determined.  Optional; default is -1 
- <strong>`Layers`</strong>:  .  Optional 
    - <strong>`Patch Name`</strong>:  Name of the patch, defined in the `Geometry Definition` section.  
    - <strong>`nSurfaceLayers`</strong>:  Number of layers per patch.  Optional; default is 1 
    - <strong>`expansionRatio`</strong>:  Expansion ratio.  Optional; default is 1. 
    - <strong>`finalLayerThickness`</strong>:  Thickness (relative or absolute) of the layer furthest from the wall.  Optional; default is 1. 
    - <strong>`firstLyrThickness`</strong>:  Thickness (relative or absolute) of the layer closest to the wall.  Optional; default is -1. 
    - <strong>`thickness`</strong>:  Overall thickness of the cells.  Optional; default is -1 
    - <strong>`minThickness`</strong>:  Minimum overall thickness of the layers.  Optional; default is 1 

Only certain combinations of layer options can be used:
- `"expansionRatio"` and `"finalLayerThickness"`
- `"expansionRatio"` and `"firstLayerThickness"`
- `"thickness"` and `"finalLayerThickness"`
- `"thickness"` and `"firstLayerThickness"`
- `"expansionRatio"` and `"thickness"`

\subsection mesh_quality Mesh Quality Controls

- <strong>`maxNonOrtho`</strong>:  Specifies the maximum allowable angle made by the vector between adjacent cell centers across the common face and face normal.  Optional; default is 65 
- <strong>`maxBoundarySkewness`</strong>:  Sets the maximum allowable skewness in the boundaries.  Optional; default is 20. 
- <strong>`maxInternalSkewness`</strong>:  Sets the maximum allowable internal skewness.  Optional; default is 4. 
- <strong>`maxConcave`</strong>:  Sets the maximum concavity.  Optional; default is 80. 
- <strong>`minVol`</strong>:  Sets the minimum cell volume in cubic meters.  Optional; default is 1e-13 
- <strong>`minTetQuality`</strong>:  Sets the minimum tetrahedral element quality. Set to a small positive number (1e-15), it enables tracking. Set to a large negative number (-1e30) it enables optimal layer insertion.  Optional; default is 1e-15 
- <strong>`minArea`</strong>:  Sets the minimum area for a cell face, in square meters.  Optional; default is -1. (no minimum) 
- <strong>`minTwist`</strong>:  Sets the minimum value for twist.  Optional; default is 0.02 
- <strong>`minFaceWeight`</strong>:  Sets the minimum face interpolation weight.  Optional; default is 0.05 
- <strong>`minVolRatio`</strong>:  Sets the minimum volume ratio.  Optional; default is 0.01 
- <strong>`minDeterminant`</strong>:  Sets the minimum cell determinant.  Optional; default is 0.001 
- <strong>`minTriangleTwist`</strong>:  Sets the minimum triangle twist.  Optional; default is -1. (no minimum) 
- <strong>`qcnSmoothScale`</strong>:  Smoothing iterations, used in conjunction with `errorReduction`.  Optional; default is 5 
- <strong>`errorReduction`</strong>:  Error reduction, used in conjunction with `qcnSmoothScale`. Sets the amount to scale back displacement at points of error.  Optional; default is 0.75 
