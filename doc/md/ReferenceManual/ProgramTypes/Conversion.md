@page conversion_ref Mesh Conversion

<em>Promesh</em> provides conversion between mesh formats.

A conversion JSON file contains two sections: `"Mesh File Options"`, where the names of the input and output meshes are provided, and `"Conversion Options"` where the conversion method is specified. Most conversion methods also have their own additional requirements.

    {
        "Program Type": "Conversion",
        "Mesh File Options": {
            "Input Mesh File": "input.msh",
            "Output Mesh File": "output.vtu",
        },
        "Conversion Options": {
            "Method": "GMSH->VTK"
        }
    }

Mesh conversion options include:
- \ref gmsh_exo `"GMSH->EXO"`
- \ref gmsh_vtk `"GMSH->VTK"`
- \ref vtk_cobalt `"VTK->COBALT"`
- \ref vtk_patran `"VTK->PATRAN"`
- \ref vtk_pnt `"VTK->PNT"` 
- \ref vtk_hex_vtk_tet `"VTK-HEX->VTK-TET"` 

If CFMesh is available, further options are:
- \ref oF_msh `"FOAM->MSH"` 
- \ref oF_vtk `"FOAM->VTK"` 
- \ref vtk_oF `"VTK->FOAM"` 

\subsubsection gmsh_exo gmsh to EXODUS II

Within the `"Conversion Options"`, requires:
- <strong>`Number of Mesh`</strong>: the number of meshes (must be at least one)
- <strong>`Mesh Data`</strong>: block which contains:
    - <strong>`File`</strong>: the mesh's file name <em>(Required)</em>
    - <strong>`Name`</strong>: the mesh's name <em>(Optional, default `"default"`)</em>
    - <strong>`Use Physical Groups`</strong>: can be `true` or `false` <em>(Optional, default `false`)</em>
    - <strong>`Node Sets`</strong>: gives an array listing the node sets, if any
    - <strong>`Element Blocks`</strong>: gives an array listing the element blocks, if any
    - <strong>`Free Surface Side Set`</strong>: can be `true` or `false` <em>(Optional, default `false`)</em>
    - <strong>`Split Top and Bottom`</strong>: can be `true` or `false` <em>(Optional, default `false`)</em>
    - <strong>`Side Set Names`</strong>: <em>(Optional)</em> 
    - <strong>`Element Block Names`</strong>: <em>(Optional)</em>
    - <strong>`Add Global Node Set`</strong>: <em>(Optional)</em>
- <strong>`Post Processing`</strong>: requires a post processing file in JSON format to be specified
- <strong>`Number of Tasks`</strong>: <em>(default 0)</em>
- <strong>`Tasks`</strong>: blocks should correspond with the number of tasks input above
    - <strong>`File`</strong>:gives the name of the task file.

An example `"GMSH->EXO"` JSON file named `gmsh2exo.json` can be found in the `testing/test_data/ep16pre` directory.

\subsubsection gmsh_vtk gmsh to VTK
Requires no additional keywords or specifications.

\subsubsection vtk_cobalt VTK to COBALT
Within the `"Mesh File Options"`, requires:
- <strong>`Output Patch Map File`</strong>: which is a `.cgi` file.

\subsubsection vtk_patran VTK to PATRAN

Within the `"Conversion Options"`, requires:
- <strong>`BC Info`</strong>:
    - <strong>`Patch Number`</strong>: corresponds to numbers assigned to patched within the VTK file
    - <strong>`BC Type`</strong>: can be either: `"Face"` or `"Node"`, which changes the following options as described:
    - <strong>`RocFrac FSI Type`</strong>: can be 0, 1, or 2  <em>(for `"BC Type": "Face"`)</em>
    - <strong>`Structural`</strong>: can be `true` or `false` <em>(for `"BC Type": "Node"`)</em>
    - <strong>`Mesh Motion`</strong>: can be `true` or `false` <em>(for `"BC Type": "Node"`)</em>
    - <strong>`Thermal`</strong>: can be `true` or `false` <em>(for `"BC Type": "Node"`)</em>
    - <strong>`RocfracControl Type`</strong>: can be 1, 2, or 3 <em>(for `"BC Type": "Node"`)</em>
- <strong>`Node Patch Preference`</strong>: lists the preferred order for the node patches

An example `"VTK->PATRAN"` JSON file named `vtk2patran.json` can be found in the `testing/test_data/PatranGenTest` directory.

\subsubsection vtk_pnt VTK to PNT
Within the `"Conversion Options"`, requires:
- <strong>`Check Conversion Quality`</strong>: can be either `true` or `false`
- <strong>`Dimension`</strong>: can be either 1, 2, or 3
- <strong>`Block Data`</strong>: a JSON file may have one or more of these sections, depending on the number of blocks present in the mesh
    - <strong>`Element Order`</strong>: can be 1 or 2
    - <strong>`Equation Order`</strong>: can be 1 or 2
    - <strong>`Element Type`</strong>: can be `"BAR"` for 1D, `"TRIANGLE"`, `"QUADRILATERAL"`, or `"HEXAGON"` for 2D, or `"SPHERICAL"`, `"CYLINDRICAL"`, `"BRICK"`, `"LAGRANGE_BRICK"`, `"TETRAHEDRON"`, `"HEXPRISM"`, or `"PRISMATIC"` for 3D. `"OTHER"` is also an option.
    - <strong>`BC Tag`</strong>: can be `"VOID"` `"REFLECTIVE"`
    - <strong>`Name`</strong>: is the name of the region
    - <strong>`Element ID Range`</strong>: can be two integers enclosed in square brackets, ie [0, 5], representing all integers from 0 to 5
    - <strong>`Element ID`</strong>: can be an array of integers enclosed in square brackets, ie [0, 1, 2]

Note that either `"Element ID Range"` or `"Element ID"` should be specified for each block.

Example `"VTK->PNT"` JSON files named `bench1.json`, `bench5.json`, and `bench6.json` can be found in the `testing/test_data/PNTGenTest` directory.

\subsubsection vtk_hex_vtk_tet VTK Hexahedral to VTK Tetrahedral
Requires no additional keywords or specifications.

\subsubsection oF_msh openFOAM to MSH
Requires no additional keywords or specifications.

\subsubsection oF_vtk openFOAM to VTK
Requires no additional keywords or specifications.

\subsubsection vtk_oF VTK to openFOAM
Requires no additional keywords or specifications.
