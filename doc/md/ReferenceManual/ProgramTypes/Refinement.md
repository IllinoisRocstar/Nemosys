@page refinement_ref Refinement

<em>Adaptive Mesh Refinement</em> is a meshing service designed to apply 
on-the-fly mesh refinement and coarsening for CFD and finite element modeling 
simulations based on machine learning and solution field values. This service 
currently has two modules, one for hexahedral meshes and one for tetrahedral 
meshes, both of which support user-defined and machine-learning-based 
refinement and coarsening.

The image below shows an example of refinement using <em>Omega_h</em>. 
In this case, the refined mesh (on the right) has been significantly coarsened 
compared to the unrefined mesh (left), especially in the far end where stress 
(denoted with color) is lower.
![A meshed beam before (left) and after (right) refinement with Omega h](Refinement/beam_refinement.png)

<strong>Refinement JSON Template</strong>

A mesh generation JSON file contains two sections: `"Mesh File Options"`, where 
the names of the input and output meshes are specified, and `"Refinement Options"`, 
where the refinement method and its associated parameters are specified.
```json lines
{
  "Program Type": "Refinement",
  "Mesh File Options": {
    "Input Mesh File": "unrefined_mesh.vtu",
    "Output Mesh File": "refined_mesh.vtu"
  },
  "Refinement Options": {
    "Refinement Method": "<REFINEMENT METHOD>",
    ...
  }
}
```
`"<REFINEMENT METHOD>"` can be replaced with one of the refinement methods 
available within <em>Promesh</em>:
- `"Omega_h"`
- `"value"` (uses <em>OpenFOAM</em>)
- `"gradient"` (uses <em>OpenFOAM</em>)
- `"uniform"` (uses <em>Gmsh</em>)
- `"FV"`
- `"Z2 Error Estimator"` (uses <em>Gmsh</em>)


\subsection omega_h Omega_h

- <strong>`Metric Sources`</strong>:  Defines the set of metric tensors that 
<em>Omega_h</em> will intersect to form the size field that the mesh refinement 
algorithm targets.  Required 
  - <strong>`Type`</strong>:  One of `"Constant"`, `"Variation"`, 
  `"Derivative"`, `"Given"`, `"Implied"`, or `"Curvature"`  Required; see below 
  for further explanation
  - <strong>`Knob`</strong>:  Scales the size field by multiplying the metric 
  by the inverse square of `"Knob"`.  Required
  - <strong>`Tag Name`</strong>:  Required for `"Variation"`, `"Derivative"`, 
  and `"Given"` types, ignored for other types 
  - <strong>`Isotropy`</strong>:  One of `"Anisotropic"`, `"Length"`, or 
  `"Size"`.  Required. Default is `"Anisotropic"`
  - <strong>`Scales`</strong>:  Either `"Absolute"` or `"Scales"`.  Required. 
  Default is `"Scales"` 
- <strong>`Reconstruct Geometry`</strong>:    Optional; default is `false` 
- <strong>`Verbose`</strong>:    Optional 
- <strong>`Should Limit Lengths`</strong>:  If `true`, allows the size field to 
be clamped so that the edge lengths in metric space are between `"Min Length"` 
and `"Max Length"`.  Optional 
- <strong>`Max Length`</strong>:  Sets the maximum edge length in metric space 
if `"Should Limit Lengths"` is `true`.  Optional 
- <strong>`Min Length`</strong>:  Sets the minimum edge length in metric space 
if `"Should Limit Lengths"` is `true`.  Optional 
- <strong>`Should Limit Gradation`</strong>:  If `true`, adjusts the metric 
size field to limit its gradation.  Optional 
- <strong>`Max Gradation Rate`</strong>:  Sets the maximum gradation if 
`"Max Gradation Rate"` is `true`.  Optional 
- <strong>`Gradation Convergence Tolerance`</strong>:  Sets the gradation 
convergence tolerance if `"Max Gradation Rate"` is `true`.  Optional 
- <strong>`Should Limit Element Count`</strong>:  If `true`, allows the size 
field to be scaled so that the resulting number of elements falls between 
`"Min Element Count"` and `"Max Element Count"`.  Optional 
- <strong>`Max ElementCount`</strong>:  Sets the maximum element count if 
`"Should Limit Element Count"` is `true`.  Optional 
- <strong>`Min Element Count`</strong>:  Sets the minimum element count if 
`"Should Limit Element Count"` is `true`.   Optional 
- <strong>`Element Count Over Relaxation`</strong>:  Describes how aggressively
the metric should be scaled.  Optional 
- <strong>`N smoothing Steps`</strong>:  Specifies the number of times the 
metric smoothing algorithm should be applied.  Optional 
- <strong>`Min Length Desired`</strong>:  Minimum length (in metric space) 
before edge is a candidate for coarsening.  Optional; default is \f$1/\sqrt{2}\f$ 
- <strong>`Max Length Desired`</strong>:  Maximum length (in metric space) 
before edge is a candidate for refinement.  Optional; default is \f$\sqrt{2}\f$ 
- <strong>`Max Length Allowed`</strong>:  Maximum length (in metric space) 
above which the algorithm will continue to iterate.  Optional; default is 
\f$2/\sqrt{2}\f$ 
- <strong>`Min Quality Allowed`</strong>:  Minimum quality below which the 
algorithm will continue to iterate.  Optional; default is 0.3 for 2D meshes and 
0.2 for 3D meshes 
- <strong>`Min Quality Desired`</strong>:  Specifies the threshold below which 
elements are marked as slivers for coarsening/edge swapping  Optional; default 
is 0.4 for 2D meshes and 0.3 for 3D meshes 
- <strong>`N sliver Layers`</strong>:    Optional 
- <strong>`Verbosity`</strong>:    Optional 
- <strong>`Length Histogram Min`</strong>:    Optional 
- <strong>`Length Histogram Max`</strong>:    Optional 
- <strong>`N length Histogram Bins`</strong>:    Optional 
- <strong>`N quality Histogram Bins`</strong>:    Optional 
- <strong>`Should Refine`</strong>:    Optional 
- <strong>`Should Coarsen`</strong>:    Optional 
- <strong>`Should Swap`</strong>:    Optional 
- <strong>`Should Coarsen Slivers`</strong>:    Optional 
- <strong>`Should Prevent Coarsen Flip`</strong>:    Optional 
- <strong>`Transfer Options`</strong>:    Optional 
- <strong>`Name`</strong>:  Specifies the name of a data array/tag.  Required 
  - <strong>`Method`</strong>:  One of `"Inherit"`, `"Linear Interp"`, 
  `"Metric"`, `"Pointwise"`, `"Density"`, or `"Conserve"`.  Required. Only 
  `"Linear Interp"`, `"Inherit"`, and `"Pointwise"` are meaningful for 
  quad/hex meshes 
  - <strong>`Integral Name`</strong>:    Required if transfer method is 
  `"Conserve"` 
  - <strong>`Transfer Integral Options`</strong>:    Required if transfer 
  method is `"Conserve"` 
      - <strong>`Integral Name`</strong>:    Required 
      - <strong>`Type`</strong>:    Required 
      - <strong>`Tolerance`</strong>:    Optional 
      - <strong>`Floor`</strong>:    Optional 

<strong>Metric Sources Options</strong>

| Type                          | Description                                                                                                                   |
|-------------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| <strong>`Constant`</strong>   | Specifies a uniform size field                                                                                                |
| <strong>`Variation`</strong>  | Specifies a size field based on the derivative of `Tag Name`, which must be the name of a tag on the elements of the vertices |
| <strong>`Derivative`</strong> | Specifies a size field based on the Hessian of `Tag Name`, which must be the name of a tag on the elements of the vertices    |
| <strong>`Given`</strong>      | Specifies a size field given by `Tag Name`, which must be the name of a tag on the vertices                                   |
| <strong>`Implied`</strong>    | Specifies the size field implied by the current mesh                                                                          |
| <strong>`Curvature`</strong>  | Specifies a size field based on the curvature of the surface of the mesh                                                      |



| Isotropy                        | Description                                                                                                                   |
|---------------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| <strong>`Anisotropic`</strong>  | Does not alter the metric source                                                                                              |
| <strong>`Length`</strong>                          | Turns an anisotropic metric into an isotropic one by choosing the maximum metric in any direction                             |
| <strong>`Size`</strong>                           | Turns an anisotropic mesh into an isotropic one by choosing the isotropic mesh that preserves the area/volume of the elements |

| Scales | Description                                                                              |
|--------|------------------------------------------------------------------------------------------|
| <strong>`Absolute`</strong>       | Prevents scaling of the metric source to compensate for a specific targeted element count |
| <strong>`Scales`</strong>       |  Allows scaling of the metric source to compensate for a specific targeted element count|


<strong>Transfer Options</strong>

| Method | Description                                                                                                                                       |
|--------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| <strong>`Inherit`</strong>       | Copies values from higher dimensional elements to lower dimensional ones. It is assumed that there exists a tag named `"Name"` for all dimensions |
| <strong>`Linear Interp`</strong>       | Linearly interpolates `"Name"`, assumed to be defined on points, to the new mesh                                                                  |
| <strong>`Metric`</strong>       | Interpolates `"Name"`, defined on the vertices, by treating it as a metric tensor                                                                 |
| <strong>`Pointwise`</strong>       | Copies value from the nearest cell on the original mesh                                                                                           |
| <strong>`Density`</strong>       | Conserves some quantity when transferring                                                                                                         |
| <strong>`Conserve`</strong>       | Conserves some quantity when transferring. Requires specifying `"Integal Name"`                                                                   |

  
<strong>Transfer Integral Options</strong>

| Type | Description |
|------|-------------|
| <strong>`None`</strong>     |             |
| <strong>`Relative`</strong>     |             |
| <strong>`Absolute`</strong>     |             |


\subsection size_field Size Field Refinement

For <em>Gmsh</em> based refinement.

This is indicated with:
```json lines
  "Refinement Options": {
    "Refinement Method": "value",
    ...
  }
```
or
```json lines
  "Refinement Options": {
    "Refinement Method": "gradient",
    ...
  }
```

but the subsequent parameters are the same.
- <strong>`Array Name`</strong>:    Required 
- <strong>`StdDev Multiplier`</strong>:    Required 
- <strong>`Max Is Min for Scaling`</strong>:    Required 
- <strong>`Transfer Data`</strong>:    Required 
- <strong>`Size Factor`</strong>:    Optional; default is 1. 


\subsection unif Uniform

- <strong>`Edge Scaling`</strong>:     
- <strong>`Transfer Data`</strong>:     


\subsection fv FV

<em>OpenFOAM</em> machine learning based refinement.

```json lines
{
  "Program Type": "Refinement",
  "Mesh File Options": {
    "Input Mesh File": "unrefined_mesh.vtu",
    "Output Mesh File": "refined_mesh.vtu"
  },
  "Refinement Options": {
    "Refinement Method": "<REFINEMENT METHOD>",
    ...
  }
}
```
- <strong>`Input Field File`</strong>: Provides the file containing the scalar field of the refinement criteria quantity. Required.
- <strong>`Refinement Interval`</strong>: Refinement interval/frequency. Required; default is 1. 
- <strong>`Maximum Refinement`</strong>:  Sets the maximum number of times a cell can be refined. Required; default is 1. 
- <strong>`Buffer Layers`</strong>: Buffer layers for slower than 2:1 refinement. Optional; default is 1.
- <strong>`Max Cells`</strong>: Sets the maximum number of cells; if this limit is reached, refinement will stop. Optional; default is 500000.
- <strong>`Write Field Data?`</strong>:    default is `false` 
- <strong>`Write Mesh?`</strong>:    default is `false` 
- <strong>`Write Refinement Data?`</strong>:    default is `false` 
- <strong>`Time Step`</strong>:     
- <strong>`End Time`</strong>:     
- <strong>`Start Time`</strong>:  
- <strong>`Lower Refinement Level`</strong>: Sets the lower refinement level for the refinement criteria quantity. Optional.
- <strong>`Upper Refinement Level`</strong>: Sets the upper refinement level for the refinement criteria quantity. Optional.
- <strong>`Unrefinement Above`</strong>: Sets the criteria for unrefinement/coarsening above the specified value. Optional.
- <strong>`Unrefinement Below`</strong>: Sets the criteria for unrefinement/coarsening below the specified value. Optional.


- <strong>`ML Kernel`</strong>:     

\subsection z2ee Z2 Error Estimator

- <strong>`ImproveMeshQuality`</strong>:    Required 
- <strong>`Transfer Data`</strong>:    Required 
- <strong>`Array Name`</strong>:    Required 
- <strong>`Order`</strong>:    Required 