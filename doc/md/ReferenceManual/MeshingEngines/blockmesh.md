@page blockmesh_ref blockMesh

[TOC]

\section desc Description

<em>blockMesh</em> is a meshing engine used to create hex meshes in elementary geometrical shapes: 
- \ref rprism 
- \ref spheres
- \ref cyl_cone

It is included with <em>OpenFOAM</em> and requires that the environment has been properly set. On a default
installation, this can be done by executing the following command in the working terminal:

    . /usr/lib/openfoam/openfoam2006/etc/bashrc

It is also possible to create assemblies of these elementary shapes using an input dictionary via the `"Input Dict File"` keyword. For more information, see the [<em>OpenFOAM</em> User Guide](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.3-mesh-generation-with-the-blockmesh-utility).


\section ex Examples


Generic code for <em>blockMesh</em>:

    {
        "Program Type": "Mesh Generation",
        "Mesh File Options": {
            "Output Mesh File": "filename.vtu"
        },
        "Mesh Generation Options": {
            "Mesh Generation Engine": "blockMesh",
            "blockMesh Parameters": {
                "Input Dict File": false,
                "scaleToMeters": 1,
                "NumCells": [40, 40, 40],
                "<SHAPE> Parameters": {
                    ...
                }
            }
        }
    }


where `"<SHAPE>"` can be replaced with `"Block"`, `"Sphere"`, or `"Cylinder/Tapered_Cone"`

\section params Parameters

\subsection general_params General

- <strong>`Input Dict File`</strong>:  Specifies the input file, if any  `false` indicates no input dict should be used. 
- <strong>`scaleToMeters`</strong>:  Scales all input lengths from the default of meters.To set to millimeters, set to 0.001
- <strong>`NumCells`</strong>:  Number of cells in each direction [<em>x</em>, <em>y</em>, <em>z</em>].For 40 cells in every direction, set to [40, 40, 40]


Optionally, `"Cell_Size"` can be provided within the `"blockMesh Parameters"` block. Then, the number of cells in each direction will be determined with
\f$\frac{side length}{cell size}\f$. Note that this will override any specified values for `"NumCells"`.

\subsection shape_params Shape Specification


The `"blockMesh Parameters"` block should contain <em>exactly one</em> of `"Block Parameters"`, `"Sphere Parameters"`, or `"Cylinder/Tapered_Cone Parameters"`. 

----


\subsubsection rprism Rectangular Prisms 


![An example rectangular prism created and meshed with blockMesh](blockMsh_Box.png "An example rectangular prism created and meshed with blockMesh") 

    "Block Parameters": {
        "InitialPoint": [0., 0., 0.],
        "Length": [1., 1., 1.],
        "Grading": [1., 1., 1.]
    }


- <strong>`InitialPoint`</strong>:  Location of reference point, [<em>x</em>, <em>y</em>, <em>z</em>]For an initial point at the origin, set to [0., 0., 0.]
- <strong>`Length`</strong>:  Side length in each direction [<em>Lx</em>, <em>Ly</em>, <em>Lz</em>].For a 1x1x1 cube, set to [1., 1., 1.]
- <strong>`Grading`</strong>:  Ratio between the widths of the ending cell and the starting cell in each direction [<em>Gx</em>, <em>Gy</em>, <em>Gz</em>].For uniform cells, set to [1., 1., 1.]

----

\subsubsection spheres Spheres


![An example sphere created and meshed with blockMesh](blockMsh_Sphere.png "An example sphere created and meshed with blockMesh")

    "Sphere Parameters": {
        "Center": [1., 1., 1.],
        "Radius": 1.,
        "Grading": [1.5, 1.5, 1.5]
    }


- <strong>`Center`</strong>:  Location of center point, [<em>x</em>, <em>y</em>, <em>z</em>]  
- <strong>`Radius`</strong>:  Radius of the sphere.  
- <strong>`Grading`</strong>:  Ratio between the widths of the ending cell and the starting cell in each direction [<em>Gx</em>, <em>Gy</em>, <em>Gz</em>].   

----

\subsubsection cyl_cone Cylinders/Truncated Cones


\note This is written as `"Cylinder/Tapered_Cone"` in the source code, and should be specified as shown in the input file. Nevertheless, what is <em>meant</em> is truncated, not tapered. 


    "Cylinder/Tapered_Cone Parameters": {
        "Center": [1., 1., 1.],
        "Radius1": 1.5,
        "Radius2": 0.6,
        "Height": 2.5
        "Grading": [1., 1., 1.]
    }


- <strong>`Center`</strong>:  Location of center point, [<em>x</em>, <em>y</em>, <em>z</em>] This refers to the center of the base circle. 
- <strong>`Radius1`</strong>:  The bottom radius of the cylinder/truncated cone.  For a cylinder, the two radii should be the same.
- <strong>`Radius2`</strong>:  The top radius of the cylinder/truncated cone. For a cylinder, the two radii should be the same.
- <strong>`Height`</strong>:  The height of the cylinder/truncated cone.   
- <strong>`Grading`</strong>:  Ratio between the widths of the ending cell and the starting cell in each direction [<em>Gx</em>, <em>Gy</em>, <em>Gz</em>].For uniform cells, set to [1., 1., 1.]


----

\section see_also See Also


See the [header file](blockMeshParams_8H_source.html) for source.

