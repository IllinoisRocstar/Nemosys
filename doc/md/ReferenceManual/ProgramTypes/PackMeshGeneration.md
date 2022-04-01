@page packmesh_ref Pack Mesh Generation

[TOC]

<em>Pack Mesh</em> generates meshes out of <em>Rocpack</em> geometry files.
Currently, <em>Pack Mesh</em> is only available in GPL versions of 
<em>Promesh</em>, since it requires <em>OpenFOAM</em> for the 
hexahedral volume meshing and <em>Gmsh</em> for the surface meshing. To enable
<em>Pack Mesh</em>, you must set the CMake options `ENABLE_GPL`, 
`ENABLE_CFMESH`, and `ENABLE_GMSH` to `ON` during the build step in the
installation process.

> <strong>Note</strong>: when running <em>Promesh</em> modules that rely
> on <em>OpenFOAM</em>, the environment must be set. On a default
> installation, this can be done by executing the following command in the
> working terminal:
> 
> `. /usr/lib/openfoam/openfoam2006/etc/bashrc`

There are three workflows embedded within the <em>Pack Mesh</em> engine to 
generate 3D hex(-dominant) multi-material meshes, 3D periodic tetrahedral 
multi-material meshes, and 2D surface meshes using input pack geometries in the 
form of a <em>Rocpack</em> geometry file. The basic workflow is as follows:
1. Generate a background volume mesh.
2. Snap the extracted surfaces onto the background mesh to create conformal interfaces.
3. Split the regions into different `cellZones`.
4. Merge all the different volumes into two foam meshes: the packs and the surrounding medium.
5. Create patches for the multiple pack regions.
6. Convert the foam mesh to VTK.
7. Write mesh quality statistics to a text file.


<em>Pack Mesh</em> uses <em>snappyHexMesh</em> and <em>blockMesh</em> to 
generate meshes of packs in surrounding material, and <em>Gmsh</em> to
generate surface meshes.

## Hexahedral Pack Mesh Options

For <em>snappyHexMesh</em> and <em>blockMesh</em>, only parameters that are
required or are specific to <em>Pack Mesh</em> are included here. For other
options, check the reference manual pages for those meshing engines.

<strong>Hex Pack Mesh JSON Template</strong>

```json lines
{
  "Program Type": "Pack Mesh Generation",
  "Mesh File Options": {
    "Input Rocpack File":"packs_roc",
    "Output Pack Mesh File": "geom_pack_mesh.vtu",
    "Output Surrounding Mesh File": "geom_surrounding_mesh.vtu",
    "Output Combined Mesh File": "packmesh.vtu"
  },
  "Pack Mesh Options": {
    "Type": "hexahedral",
    "Engine": "packmesh",
    "snappyHexMesh Parameters": {
      ...
    },
    "blockMesh Parameters": {
      ...
    },
    "MeshManipulation mergeMeshes Parameters": {
      ...
    },
    "MeshManipulation createPatch Parameters": {
      ...
    }
  }
}
```

### Hex Mesh Files

Four files are referenced in the `"Mesh File Options"` section of the pack mesh generation JSON input file:
 - <strong>`Input Rocpack File`</strong>: the <em>Rocpack</em> file from which the pack mesh will be generated.
 - <strong>`Output Pack Mesh File`</strong>: the mesh of the packed objects.
 - <strong>`Output Surrounding Mesh File`</strong>: the mesh of the surrounding medium.
 - <strong>`Output Combined Mesh File`</strong>: the combined pack and surrounding medium meshes.
 
### Hex Pack Mesh Options

- <strong>`Type`</strong>:  Should be set to `"hexahedral"`.  Required 
- <strong>`Engine`</strong>:  `"packmesh"`  Required 
- <strong>`snappyHexMesh Parameters`</strong>:    Required 
    - <strong>`Castellated Mesh`</strong>:  If `true`, enables castellated mesh option.  Required 
    - <strong>`Snapping`</strong>:  If `true`, enables surface snapping option. Required 
    - <strong>`Layer Addition`</strong>:  If `true`, enables adding layers. Required 
    - <strong>`Geometry Definition`</strong>:    Required 
 `Enable Multi Patches`</strong>:  Should be enabled if the geometry file has more than one patch.  Required 
 `Input Patch Name`</strong>:  If specified, will be assigned to the entire surface of the provided geometry.  Required 
    - <strong>`Castellated Mesh Controls`</strong>:    Required 
 `CellZones`</strong>:  If `true`, indicates that there will be multiple disconnected regions within the mesh.  Required 
 `RegionRefine`</strong>:  If `true`, indicates additional refinement regions.  Required 
 `SurfaceRefine`</strong>:  If `true`, indicates additional refinement surfaces.  Required 
- <strong>`blockMesh Parameters`</strong>:    Required 
    - <strong>`Input Dict File`</strong>:  Specifies the input file, if any. `false` indicates no input dict should be used. Required 
    - <strong>`scaleToMeters`</strong>:  Scales all input lengths from the default of meters.To set to millimeters, set to 0.001
    - <strong>`NumCells`</strong>:  Number of cells in each direction [<em>x</em>, <em>y</em>, <em>z</em>].For 40 cells in every direction, set to [40, 40, 40]
    - <strong>`Block Parameters`</strong>:   Required
 `InitialPoint`</strong>:  Location of reference point, [<em>x</em>, <em>y</em>, <em>z</em>] Required if not setting `"Auto_Generate"` 
 `Length`</strong>:  Side length in each direction [<em>Lx</em>, <em>Ly</em>, <em>Lz</em>]. Required if not setting `"Auto_Generate"` 
 `Auto_Generate`</strong>:  Automatically generate the bounding box.  Required if `"InitialPoint"` and `"Length"` are not specified 
     `Offset`</strong>:   Required if `"InitialPoint"` and `"Length"` are not specified 
- <strong>`MeshManipulation mergeMeshes Parameters`</strong>:  Options for combining the pack geometries to bring them under one mesh.  Optional 
    - <strong>`Master Region Path`</strong>:  Specifies the input file, if any. `false` indicates no input dict should be used. Optional 
    - <strong>`Add Region Path`</strong>:  Specifies the input file, if any. `false` indicates no input dict should be used. Optional 
    - <strong>`overwrite?`</strong>:  Should always be set to `true` for multi-material meshing.  Optional 
- <strong>`MeshManipulation createPatch Parameters`</strong>:  Sets names for the patches used.  Optional 
    - <strong>`Surrounding PatchName`</strong>:  Specifies the name for the surrounding medium. Optional 
    - <strong>`Packs PatchName`</strong>:  Specifies the name for the packed material.  Optional 
    - <strong>`Surrounding PatchType`</strong>:  Specifies the patch type (either `patch` or `wall`) for the surrounding medium.  Optional 
    - <strong>`Packs PatchType`</strong>:  Specifies the patch type (either `patch` or `wall`) for the packed material.  Optional 
    - <strong>`overwrite?`</strong>:  Should always be set to `true` for multi-material meshing.  Optional 


## Surface Pack Mesh Options

<strong>Surface Pack Mesh JSON Template</strong>
```json lines
{
  "Program Type": "Pack Mesh Generation",
  "Mesh File Options": {
    "Input Rocpack File":"packs_roc",
    "Output Mesh File": "geom_pack_mesh.vtu",
  },
  "Pack Mesh Options": {
    "Type": "surface",
    "Engine": "packmesh",
    "Periodic 3D Mesh": {
      "Physical Group Options": "None",
      "Create cohesive elements": false,
      "Enable Patches": false,
      "Set Periodic Geometry": false,
      "Element Order": 1,
      "Custom Domain": {
        "Initial": [0.0, 0.0, 0.0],
        "Length": [1.0, 1.0, 1.0]
      },
      "Transfer Mesh": [0.5, 0.5, 0.5]
    },
    "Mesh Size": 1.0,
    "Mesh Algorithm": 1,
    "Scale Value": 1,
    "Remove geometries on boundary": false,
    "Enable Default Outputs": false,
    "Enable Size Preservation": false,
    "Refinement Levels":,
    "Upper Threshold": ,
    "Lower Threshold":
  }
}
```

### Surface Mesh Files

Four files are referenced in the `"Mesh File Options"` section of the pack mesh generation JSON input file:
 - <strong>`Input Rocpack File`</strong>: the <em>Rocpack</em> file from which the pack mesh will be generated.
 - <strong>`Output Mesh File`</strong>: the mesh of the packed objects.

### Surface Pack Mesh Options

- <strong>`Type`</strong>:  Should be set to `"surface"`. Required 
- <strong>`Engine`</strong>:  `"packmesh"`  Required 
- <strong>`Periodic 3D Mesh`</strong>:    Required 
    - <strong>`Physical Group Options`</strong>:  One of `"None"`, `"Multi Physical Groups"`, `"Two Physical Groups"`, or `"Physical Group per Shape"`  Optional; default is `"None"` 
    - <strong>`Create cohesive elements`</strong>:  If `true`, creates zero-thickness 3D elements between the pack and the surrounding region.  Optional; default is `false` 
    - <strong>`Enable Patches`</strong>:  If `true`, enables writing of surface patches in the final mesh.  Optional; default is `false` 
    - <strong>`Set Periodic Geometry`</strong>:    Optional; default is `false` 
    - <strong>`Element Order`</strong>:  Either 1 or 2  Optional; default is 1 
    - <strong>`Custom Domain`</strong>:    Optional 
      - `Initial`</strong>:  Provides the [<em>x</em>, <em>y</em>, <em>z</em>] coordinates of the initial point
      - `Length`</strong>: Specifies the lengths in the [<em>x</em>, <em>y</em>, <em>z</em>] directions.
    - <strong>`TransferMesh`</strong>:  Specifies translation in the [<em>x</em>, <em>y</em>, <em>z</em>] directions.  Optional 
- <strong>`Mesh Size`</strong>:    Optional 
- <strong>`Mesh Algorithm`</strong>:  Options are 1, 2, 5, 6, 7, 8, or 9 for 2D meshes and 1, 4, 7, 9, or 10 for 3D meshes; see below.  Optional; default is 1 

|  Dimension  |  Number  | Name                       | Description                                                                                                                                                     |
|:-----------:|:--------:|----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
|     2D      |    1     | MeshAdapt                  | Uses edge swaps, splits, and collapses to obtain a better geometrical configuration                                                                             |
|      ^      |    2     | Automatic                  | Attempts to choose the best algorithm                                                                                                                           |
|      ^      |    5     | Delaunay                   | Inserts new points sequentially at the circumceter of the element with the largest circumradius, then reconnects the mesh with an anisotropic Delaunay criterion |
|      ^      |    6     | Frontal Delaunay           | Computes point positions and connections simultaneously                                                                                                         |
|      ^      |    7     | BAMG                       | Allows generating anisotropic triangulations                                                                                                                    |
|      ^      |    8     | Frontal Delaunay for Quads | Variation of Frontal Delaunay that generates right-angle triangles suitable for recombination                                                                   |
|      ^      |    9     | Packing of Parallelograms  | Fills the mesh with parallelograms that then get converted to triangles                                                                                         |
|     3D      |    1     | Delaunay                   | First creates an initial mesh out of the union of all volumes in the model, then a 3D version og the 2D Delaunay algorithm is applied                           |
|      ^      |    4     | Frontal                    | Uses the Netgen algorithm                                                                                                                                       |
|      ^      |    7     | MMG3D                      | Generates anisotropic tetrahedralizations                                                                                                                       |
|      ^      |    9     | R-Tree                     | Allows for further recombination into hexes                                                                                                                     |
|      ^      |    10    | HXT                        | Efficient, parallel reimplementation of the Delaunay algorithm                                                                                                  |

- <strong>`Scale Value`</strong>:  Sets the amount of reduction, from 0 to 1, with 1 being no scaling  Optional 
- <strong>`Remove geometries on boundary`</strong>:  If `true`, any geometry intersecting the boundary will be removed  Optional; default is `false` 
- <strong>`Enable Default Outputs`</strong>:    Optional; default is `false` 
- <strong>`Enable Size Preservation`</strong>:  If `true`, pack size is preserved instead of packing fraction.  Optional; default is `false` 
- <strong>`Refinement Levels`</strong>:  Refinement applied to the original mesh.  Optional 
- <strong>`Upper Threshold`</strong>:  Upper threshold for filtering with respect to mean.  Optional 
- <strong>`Lower Threshold`</strong>  Lower threshold for filtering with respect to mean.  Optional 
