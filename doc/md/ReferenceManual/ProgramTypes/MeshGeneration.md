@page meshgeneration_ref Mesh Generation

Meshes can be generated from existing CAD models or simple shapes specified within the input files.

<strong>Mesh Generation JSON Template</strong>

A mesh generation JSON file contains two sections: `"Mesh File Options"`, where the names of the input geometry and output mesh are specified, and `"Mesh Generation Options"`, where the meshing engine and its associated parameters are specified.

    {
        "Program Type": "Mesh Generation",
        "Mesh File Options": {
            "Input Geometry File": "geometry.stl"
            "Output Mesh File": "meshed_geometry.vtu"
        },
        "Mesh Generation Options": {
            "Mesh Generation Engine": "netgen",
            "Netgen Parameters": {
                ...
            }
        }
    }

<strong>Meshing Engines</strong>

<em>Promesh</em> has several meshing engines available:
- @subpage blockmesh_ref
- @subpage cfmesh_ref
- @subpage netgen_ref
- @subpage snappyhexmesh_ref

These tools can be used individually or in conjunction with one another.
