@page gmsh_ref Gmsh

\section desc Description

The <em>Gmsh</em> meshing engine is an interface to the well-known, commonly 
used <em>Gmsh</em> meshing tool. This pipeline mainly supports IGES and STEP
input geometries.

A generic <em>Gmsh</em> input file is shown, with the default values for each 
parameter provided. Note that all parameters listed under `"Gmsh Parameters"` 
are optional.
```json
{
  "Program Type": "Mesh Generation",
  "Mesh File Options": {
    "Input Geometry File": "geometry.step",
    "Output Mesh File": "meshed_geometry.vtu"
  },
  "Mesh Generation Options": {
    "Mesh Generation Engine": "gmsh",
    "Gmsh Parameters": {
      "minSize": 0.01,
      "maxSize": 50.0,
      "surfaceAlgorithm": "Frontal",
      "volumeAlgorithm": "HXT",
      "extendSizeFromBoundary": true,
      "sizeFromCurvature": false,
      "minElementsPer2Pi": 6,
      "optimize": false,
      "optimizeThreshold": 0.3,
      "BackgroundField": -1,
      "elementOrder": 1,
      "subdivisionAlgorithm": 1,
      "saveAll": true,
      "fragment": false,
      "TransfiniteBlocks": [
        {
          "Volume": 1,
          "Axis": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "x": {
            "Vertices": 3,
            "Progression": 1.0
          },
          "y": {
            "Vertices": 3,
            "Bump": 0.1
          },
          "z": {
            "Vertices": 3,
            "Bump": 0.1
          }
        }
      ],
      "ColorMap": [
        {
          "Color": [0, 0, 0, 0],
          "Group": "group_name"
        }
      ],
      "SizeFields": [
        {
          "Type": "Ball",
          "ID": 2,
          "Params": {
            ...
          }
        }
      ]
    }
  }
}
```

\section params Gmsh Parameters

- <strong>`minSize`</strong>: Minimum global mesh size; this value is dimensionless. 
Default is 0.01.
- <strong>`maxSize`</strong>: Maximum global mesh size; this value is dimensionless. 
Default is 50.0.
- <strong>`surfaceAlgorithm`</strong>: Surface meshing algorithm; see table below for 
the available options and their descriptions. The default is `"Frontal"`.

| Surface Meshing Algorithm   | Description                                                                                                     |
|-----------------------------|-----------------------------------------------------------------------------------------------------------------|
| <strong>`Frontal`</strong>                   | Uses an advancing front algorithm                                                                               |
| <strong>`MeshAdapt`</strong>                 | Uses edge sweeps, splits, and collapses to adapt the mesh                                                       |
| <strong>`Delaunay`</strong>                  | Uses Delaunay triangulation                                                                                     |
| <strong>`Frontal Quad`</strong>              | An experimental algorithm that attempts to create right-angle triangles, making recombination to quads smoother |
| <strong>`Packing of Parallelograms`</strong> | Fills the mesh with parallelograms that ten get converted to triangles                                          |
| <strong>`Automatic`</strong>                 | Attempts to select the best meshing algorithm based on the geometry                                             |
> Note that the `"Frontal"`, `"Frontal Quad"`, and `"Packing of Parallelograms"`
> surface algorithms are not entirely compatible with `"SizeFields"`. If a size
> field is defined that intersects with a surface with these algorithms defined,
> the surface mesh may not reflect the mesh size defined in the size field.

- <strong>`volumeAlgorithm`</strong>: Volume meshing algorithm; see table below
for the available options and their descriptions. The default is `"HXT"`.

| Volume Meshing Algorithm | Description                       |
|--------------------------|-----------------------------------|
| <strong>`Delaunay`</strong>               | Uses Delaunay triangulation       |
| <strong>`Frontal`</strong>                | Uses an advancing front algorithm |
| <strong>`HXT`</strong>                    | A fast Delaunay type              |
> Note that the `"HXT"` volume algorithm is not entirely compatible with 
> `"SizeFields"`. If a size field is defined that intersects with a 
> volume with this algorithm specified, the volume mesh may not 
> reflect the mesh size defined in the size field. 
> 
> The `"Frontal"` volume algorithm <strong>should not be used</strong> if a size
> field is used.
- <strong>`extendSizeFromBoundary`</strong>: If enabled, will extend the mesh
size from boundaries, allowing for a smooth propagation of the mesh size
smoothly from the boundary to the bulk mesh. Default is `true`.
- <strong>`sizeFromCurvature`</strong>: If enabled, will calculate the mesh
size based on geometry curvature, resulting in higher mesh generation times. 
This option is accompanied by the`"minElementsPer2Pi"` option below. Default
is `false`.
- <strong>`minElementsPer2Pi`</strong>: If `"sizeFromCurvature"` is set to
`true`, this option provides the minumum number of elements to be used per
\f$2\pi\f$ of geometric curvature. Default is 6.
- <strong>`optimize`</strong>: If set to `true`, enables volume mesh 
optimization. Default is `false`.
- <strong>`optimizeThreshold`</strong>: If `"optimize"` is set to true,
this option defines the element quality threshold for mesh optimization.
This is a number between 0 and 1, where 1 is perfect quality and 0 is an 
illegal tetrahedron. For example, setting this to 0.5 means that all mesh
elements with quality less than and equal to 0.5 will be optimized. 
Default is 0.3.
- <strong>`BackgroundField`</strong>: This option sets the ID of the size
filed to be used and is only applicable when size fields have been defined.
- <strong>`elementOrder`</strong>: Element order. Default is 1.
- <strong>`subdivisionAlgorithm`</strong>: Specifies the subdivision 
algorithm; `0` is none, `1` is all quads, and `2` is all hexes. Default is 1.
- <strong>`saveAll`</strong>: If enabled, saves all elements. Default is `true`.
- <strong>`fragment`</strong>: If enabled, calls Boolean fragments on all 
volumes. Default is `false`.
- <strong>`TransfiniteBlocks`</strong>: Provides a list of transfinite hexahedra
  (the <em>Gmsh</em> equivalent of a structured mesh).
    - <strong>`Volume`</strong>: The integer ID or tag of the volume to be 
      treated as transfinite.
    - <strong>`Axis`</strong>: Three direction vectors that specify the axis
      of the hexahedron.
    - <strong>`x`</strong>:
      - <strong>`Vertices`</strong>: The number of vertices in the <em>x</em> 
      direction.
      - <strong>`Progression`</strong>: The axis type (may be either `Progression`
      or `Bump`) and its associated coefficient. Default is `"Progression"` with
      a coefficient of 1.0.
    - <strong>`y`</strong>: 
      - <strong>`Vertices`</strong>: The number of vertices in the <em>y</em>
        direction.
      - <strong>`Bump`</strong>: A `Bump` type axis instructs the transfinite
      algorithm to distribute the vertices with a refinement at both ends of
      the line.
    - <strong>`z`</strong>:
      - <strong>`Vertices`</strong>: The number of vertices in the <em>z</em>
        direction.
      - <strong>`Progression`</strong>: A `Progression` type axis instructs the
      transfinite algorithm to distribute the vertices following a geometric
      progression: for example, `"Progression": 2.0` will result in each
      successive vertex being placed twice as far as the previous one.
- <strong>`ColorMap`</strong>: For CAD models that cannot have assigned material
groups but can have assigned colors, the `ColorMap` option provides a way to
translate between those assigned RGB values and the desired material group name.
  - <strong>`Color`</strong>: Provides an RGBA (red, green, blue, alpha) array 
with four integer values between 0 and 255.
  - <strong>`Group`</strong>: Gives the desired material group name.
- <strong>`SizeFields`</strong>: Size fields are used to control the desired mesh
size, either at geometric features (points, edges, surfaces) or in the bulk.
  - <strong>`ID`</strong>: A unique strictly positive integer used to identify
  the field.
  - <strong>`Type`</strong>: The shape of the size field; see the table below for
  the available shapes and their associated options. Additional notes are
  provided below the table.

| Type        | Parameters                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|-------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| <strong>`Ball`</strong>      | <strong>`Radius`</strong>: radius of the sphere<br/><strong>`Thickness`</strong>: the radial distance the mesh changes from VIn to VOut<br/><strong>`VIn`</strong>, <strong>`VOut`</strong>: the mesh size inside and outside the sphere <br/><strong>`XCenter`</strong>, <strong>`YCenter`</strong>, <strong>`ZCenter`</strong>: the coordinates of the sphere center                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| <strong>`Box`</strong>       | <strong>`XMin`</strong>, <strong>`YMin`</strong>, <strong>`ZMin`</strong>: minimum coordinates of the box<br/> <strong>`XMax`</strong>, <strong>`YMax`</strong>, <strong>`ZMax`</strong>: maximum coordinates of the box                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| <strong>`Cylinder`</strong>  | <strong>`Radius`</strong>: radius of the cylinder<br/> <strong>`VIn`</strong>, <strong>`VOut`</strong>: the mesh size inside and outside the cylinder<br/><strong>`XCenter`</strong>, <strong>`YCenter`</strong>, <strong>`ZCenter`</strong>: the coordinates of the cylinder center<br/><strong>`XAxis`</strong>, <strong>`YAxis`</strong>, <strong>`ZAxis`</strong>: the components of the vector defining the cylinder axis                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| <strong>`Distance`</strong>  | <strong>`EdgesList`</strong>: a list of geometry edges/curves to use for distance computation<br/><strong>`FacesList`</strong>: a list of geometry faces/surfaces to use for distance computation<br/><strong>`NNodesByEdge`</strong>: number of nodes to discretize the surface (if `FacesList` is used)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| <strong>`Frustum`</strong>   | <strong>`X1`</strong>, <strong>`Y1`</strong>, <strong>`Z1`</strong>: the coordinates of endpoint 1<br/> <strong>`X2`</strong>, <strong>`Y2`</strong>, <strong>`Z2`</strong>: the coordinates of endpoint 2<br/><strong>`R1_inner`</strong>, <strong>`R1_outer`</strong>: inner and outer radius of frustum at endpoint 1<br/><strong>`R2_inner`</strong>, <strong>`R2_outer`</strong>: inner and outer radius of frustum at endpoint 2<br/><strong>`V1_inner`</strong>,<strong>`V1_outer`</strong>: mesh size at point 1 at the inner and outer radius<br/> <strong>`V2_inner`</strong>. <strong>`V2_outer`</strong>: mesh size at point 2 at the inner and outer radius<br/><strong>`IField`</strong>: the reference cylinder field ID, if using the frustum with a cylinder field.<br/><strong>`Thickness`</strong>: the radial distance the mesh is graded from inner to outer |
| <strong>`Max`</strong>       | <strong>`FieldsList`</strong>: a list of size field IDs. The IDs on the list must be defined prior to using the `Max` field.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| <strong>`Min`</strong>       | <strong>`FieldsList`</strong>: a list of size field IDs. The IDs on the list must be defined prior to using the `Min` field.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| <strong>`Restrict`</strong>  | <strong>`VerticesList`</strong>: a list of points/vertices<br/> <strong>`EdgesList`</strong>: a list of edges/curves<br/> <strong>`FacesList`</strong>: a list of faces/surfaces<br/><strong>`RegionsList`</strong>: a list of volumes<br/><strong>`IField`</strong>: reference field to restrict                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| <strong>`Threshold`</strong> | <strong>`DistMin`</strong>: distance from entity up to which mesh size will be `LcMin`<br/><strong>`DistMax`</strong>: distance from entity after which mesh size will be `LcMax`<br/> <strong>`LcMin`</strong>: mesh size inside `DistMin`<br/><strong>`LcMax`</strong>: mesh size outside `DistMax`<br/><strong>`IField`</strong>: reference `Distance` field                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
<strong> Notes for size fields</strong>:

> For the `Cylinder` size field, there is no `Thickness` parameter,
> and thus the mesh size will suddenly jump from `VIn` to `VOut`. To help 
> grade the mesh smoothly away from the cylinder, the `Frustum` field can be 
> used in conjunction with the `Cylinder` field; in this case, the cylinder
> must be defined first.

> The `Distance` field does not vary the mesh size itself--it
> merely creates an underlying field for the `Threshold` field to work on.

> The `Max` and `Min` fields apply the maximum and minimum mesh sizes from
> an intersection of size fields. IDs in the `FieldsList` must be defined
> prior to using these fields.

> The `Restrict` field prevents a size field from being applied to 
> geometrical points, curves, surfaces, or volumes.