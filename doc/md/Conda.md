# Install from Conda #

[![Anaconda-Server Badge](https://anaconda.org/illinoisrocstar/promesh/badges/installer/conda.svg)](https://anaconda.org/illinoisrocstar/promesh)

## Conda ##
Conda is a package and environment manager that enables installation of the *Promesh*
package, [available here,](https://anaconda.org/IllinoisRocstar/promesh) with the single command :
```commandline
conda install -c illinoisrocstar promesh 
```

Instructions for installing Conda can be [found here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

## Running Promesh ##


### Creating the Input File ###

The simplest method to run a *Promesh* program is to start with a JSON input
file. A sample one is provided below:

```json
{
 "Program Type": "NucMesh Generation",
 "Mesh File Options": {
  "Output Mesh File": "test.vtu"
 },
 "NucMesh Options": {
  "Shapes": [
   {
    "Type": "Rectangular Array",
    "Grid Distance": [1.5, 1.5],
    "Pattern": [[0, 1],[1, 0]],
    "Shapes": [
     {
      "Type": "Circles",
      "Rings": [
       {
        "Radius": 0.5,
        "Mesh": {"Type": "T"},
        "Material": "A"
       },
       {
        "Radius": 1.0,
        "Mesh": {
         "Type": "S",
         "Number of Elems": [3, 8]
        },
        "Material": "B"
       }
      ]
     },
     {
      "Type": "Circles",
      "Rings": [
       {
        "Radius": 0.5,
        "Mesh": {"Type": "Q"},
        "Material": "C"
       },
       {
        "Radius": 1.0,
        "Mesh": {
         "Type": "S",
         "Number of Elems": [3, 8]
        },
        "Material": "D"
       }
      ]
     }
    ]
   }
  ]
 }
}
```

### Executing Promesh ###
To run this program, copy and paste the input above into your favorite text editor
and save as `test_run.json`. Within a terminal, navigate to the folder containing 
the JSON and enter:
```commandline
nemosysRun test_run.json
```
This should generate the following printout in the terminal:

```console
srvBase constructed
Gmsh initialized
geoMeshBase constructed
srvBase destructed
geoMeshBase constructed
vtkGeoMesh constructed
vtkGeoMesh destructed
geoMeshBase destructed
geoMeshBase destructed
Gmsh finalized
```

### Visualizing the Result ###
A `test.vtu` file should appear in the same directory. This can be visualized
with *Paraview* or your visualization software of choice:

![The test.vtu mesh generated with NucMesh](doc/images/test_run.png)
