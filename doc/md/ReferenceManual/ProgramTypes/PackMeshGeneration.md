@page packmesh_ref Pack Mesh Generation

There are three workflows embedded within the <em>Pack Mesh</em> engine to generate 3D hex(-dominant) multi-material meshes, 3D periodic tetrahedral multi-material meshes, and 2D surface meshes using input pack geometries in the form of an stl file or a <em>Rocpack</em> geometry file. The basic workflow is as follows:
- Generate a background volume mesh.
- Snap the extracted surfaces onto the background mesh to create conformal interfaces.
- Split the regions into different `cellZones`.
- Merge all the different volumes into two foam meshes: the packs and the surrounding medium.
- Create patches for the multiple pack regions.
- Convert the foam mesh to VTK.
- Write mesh quality statistics to a text file.

<strong>Pack Mesh JSON Template</strong>

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


<em>Pack Mesh</em> uses <em>snappyHexMesh</em> and <em>blockMesh</em> to generate meshes of packs in surrounding material.

\subsection pack_mesh_params Pack Mesh Options

For <em>snappyHexMesh</em> and <em>blockMesh</em>, only parameters that are required or are specific to <em>Pack Mesh</em> are included here. For other options, check the reference manual pages for those meshing engines.

\subsection hex_pack Hexahedral

\subsubsection hex_mesh_files Mesh Files

Four files are referenced in the `"Mesh File Options"` section of the pack mesh generation JSON input file:
 - <strong>`Input Rocpack File`</strong>: the <em>Rocpack</em> file from which the pack mesh will be generated. If using an stl file, instead specify with <strong>`Input Geometry File`</strong>
 - <strong>`Output Pack Mesh File`</strong>: the mesh of the packed objects.
 - <strong>`Output Surrounding Mesh File`</strong>: the mesh of the surrounding medium.
 - <strong>`Output Combined Mesh File`</strong>: the combined pack and surrounding medium meshes.

 \subsubsection hex_params Hexahedral Pack Mesh Options


- <strong>`Type`</strong>:  One of `"surface"` or `"hexahedral"`.  Required 
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


\subsection surface_pack Surface

\subsubsection surface_mesh_files Mesh Files

Four files are referenced in the `"Mesh File Options"` section of the pack mesh generation JSON input file:
 - <strong>`Input Rocpack File`</strong>: the <em>Rocpack</em> file from which the pack mesh will be generated. If using an stl file, instead specify with <strong>`Input Geometry File`</strong>
 - <strong>`Output Pack Mesh File`</strong>: the mesh of the packed objects.
 - <strong>`Output Surrounding Mesh File`</strong>: the mesh of the surrounding medium.
 - <strong>`Output Combined Mesh File`</strong>: the combined pack and surrounding medium meshes.

 \subsubsection surface_params Surface Pack Mesh Options


- <strong>`Type`</strong>:  One of `"surface"` or `"hexahedral"`.  Required 
- <strong>`Engine`</strong>:  `"packmesh"`  Required 
- <strong>`Periodic 3D Mesh`</strong>:    Required 
    - <strong>`Physical Group Options`</strong>:  One of `"None"`, `"Multi Physical Groups"`, `"Two Physical Groups"`, or `"Physical Group per Shape"`  Optional; default is `"None"` 
    - <strong>`Create cohesive elements`</strong>:  If `true`, creates zero-thickness 3D elements between the pack and the surrounding region.  Optional; default is `false` 
    - <strong>`Enable Patches`</strong>:  If `true`, enables writing of surface patches in the final mesh.  Optional; default is `false` 
    - <strong>`Set Periodic Geometry`</strong>:    Optional; default is `false` 
    - <strong>`Element Order`</strong>:  Either 1 or 2  Optional; default is 1 
    - <strong>`Custom Domain`</strong>:    Optional 
 `Initial`</strong>:    
 `Length`</strong>:    
    - <strong>`TransferMesh`</strong>:  Specifies translation in the [<em>x</em>, <em>y</em>, <em>z</em>] directions.  Optional 
- <strong>`Mesh Size`</strong>:    Optional 
- <strong>`Mesh Algorithm`</strong>:  Options are 1, 2, 5, 6, 7, 8, or 9 for 2D meshes and 1, 4, 7, 9, or 10 for 3D meshes; see below.  Optional; default is 1 
- <strong>`Scale Value`</strong>:  Sets the amount of reduction, from 0 to 1, with 1 being no scaling  Optional 
- <strong>`Remove geometries on boundary`</strong>:  If `true`, any geometry intersecting the boundary will be removed  Optional; default is `false` 
- <strong>`Enable Default Outputs`</strong>:    Optional; default is `false` 
- <strong>`Enable Size Preservation`</strong>:  If `true`, pack size is preserved instead of packing fraction.  Optional; default is `false` 
- <strong>`Refinement Levels`</strong>:  Refinement applied to the original mesh.  Optional 
- <strong>`Upper Threshold`</strong>:  Upper threshold for filtering with respect to mean.  Optional 
- <strong>`Lower Threshold`</strong>  Lower threshold for filtering with respect to mean.  Optional 

<strong>Meshing Algorithms for 2D</strong>
- MeshAdapt (1)
- Automatic (2)
- Delaunay (5)
- Frontal Delaunay (6)
- BAMG (7)
- Frontal Delaunay for Quads (8)
- Packing of Parallelograms (9)

<strong>Meshing Algorithms for 3D</strong>
- Delaunay (1)
- Frontal (4)
- MMG3D (7)
- R-Tree (9)
- HXT (10)