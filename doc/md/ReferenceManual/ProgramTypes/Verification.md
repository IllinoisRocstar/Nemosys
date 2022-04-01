@page verification_ref Verification/Mesh Convergence

Solution verification checks for convergence by testing the solution on 
progressively finer grids. Three successively refined meshes of a given 
geometry are required to evaluate whether the coarsest mesh is within the 
asymptotic range. The grid convergence index (GCI) of each component of each 
selected field is evaluated: if they are all approximately equal to 1 then the 
coarse mesh is sufficiently refined for an accurate solution (and is within 
asymptotic range).

A verification JSON file contains two sections: `"Mesh File Options"`, where 
the names of the three meshes to be used are provided, and 
`"Verification Options"`, where the solution array IDs, transfer type, and 
number of threads are provided.
```json
{
  "Program Type": "Verification",
  "Mesh File Options": {
    "Coarse Mesh File": "coarse.vtu",
    "Fine Mesh File": "fine.vtu",
    "Finer Mesh File": "finer.vtu"
  },
  "Verification Options": {
    "Array IDs": [0],
    "Transfer Type": "Consistent Interpolation",
    "Target GCI": 1.1,
    "Threads": 2
  }
}
```
Other options for `"Transfer Type"` are `"Conservative Surface Transfer"` (if 
enabled when compiling <em>Promesh</em>; requires <em>IMPACT</em>) and 
`"Conservative Volume Transfer"` (if enabled when compiling <em>Promesh</em>).  
The default GCI is 1.1, but a target GCI can be provided with keyword `"Target GCI"`. 
`"Threads"` sets the number of threads used in the transfer (if OpenMP is enabled); 
if the keyword is unset and OpenMP is enabled, @c omp_get_max_threads will be used.