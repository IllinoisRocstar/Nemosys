@page meshquality_ref Mesh Quality


Currently, <em>Promesh</em> uses <em>cfMesh</em> as a mesh quality optimizer. Because it is based on <em>openFOAM</em>, it assumes the folders that contain the mesh is present in the current directory. Therefore, there are no `"Mesh File Options"` and no mesh input or output files specified in JSON input files for `"Mesh Quality"` programs.

<strong>Mesh Quality JSON Template</strong>

A mesh quality JSON file contains a single section: `"Mesh Quality Options"`, where the meshing engine and its associated parameters are specified.

    {
        "Program Type": "Mesh Quality",
        "Mesh Quality Options": {
            "Mesh Quality Engine": "cfmesh",
            "Schedule": [
                {
                    "Method": "meshOptimizer",
                    "NIterations": 4,
                    "NLoops": 4,
                    "QualityThreshold": 0.3,
                    "NSurfaceIterations": 4,
                    "ConstrainedCellSet": null
                }
            ]
        }
    }


Mesh quality is improved via smoothing.

- <strong>`Method`</strong>:  Should be set to `"meshOptimizer"`   
    - <strong>`NIterations`</strong>:  Number of optimization iterations.  Default is 50 
    - <strong>`NLoops`</strong>:  Number of inner loops in optimization.  Default is 10 
    - <strong>`QualityThreshold`</strong>:  Minimum mesh quality before optimization is considered done.  Default is 0.1 
    - <strong>`NSurfaceIterations`</strong>:  Number of surface iterations.  Default is 2 
    - <strong>`ConstrainedCellsSet`</strong>:  Name of constrained cell set.  