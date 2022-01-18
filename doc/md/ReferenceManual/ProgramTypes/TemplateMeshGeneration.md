@page template_ref Template Mesh Generation


<strong>Template Mesh JSON Template</strong>

    {
        "Program Type": "Template Mesh Generation",
        "Mesh File Options": {
            "Output Mesh File": "spiral"
        },
        "Template Options": {
            "Template Name": "Spiral Tape Pipe",
            "Params": {
                ...
            }
        }
    }

\subsection template_params Template Parameters


- <strong>`rx`</strong>:  Ellipse radius (<em>x</em>-direction). Default is 1. 
- <strong>`ry`</strong>:  Ellipse radius (<em>y</em>-direction). Default is 1. 
- <strong>`thickness`</strong>:  Thickness of spiral tape in <em>y</em>-direction. Default is 0.25 
- <strong>`extrude_len`</strong>:  Extrusion length. Default is 3. 
- <strong>`n_turns`</strong>:  Number of spiral turns. Default is 0.5 
- <strong>`width_percent`</strong>:  Percentage of <em>x</em>-direction diameter for tape width. Default is 0.85 
- <strong>`mSize_min`</strong>:  Minimum mesh size, size at walls. Default is 0.048 
- <strong>`mSize_max`</strong>:  Maximum mesh size, size in bulk. Default is 0.1 
- <strong>`dist_min`</strong>:  Distance from walls with mSize_min default is 0.05 
- <strong>`dist_max`</strong>:  Distance from walls with mSize_max. Default is 0.2 
- <strong>`bl_wall_n`</strong>:  Boundary layer mesh size normal to wall. Default is 0.0038 
- <strong>`bl_far`</strong>:  Mesh size away from wall. Default is 0.08 
- <strong>`bl_thickness`</strong>:  Boundary layer mesh thickness. Default is 0.02 
- <strong>`ratio`</strong>:  Mesh size ratio normal to wall. Default is 1.3 
- <strong>`fan_points`</strong>:  Number of fan points in boundary layer at corners.
- <strong>`extrude_layers`</strong>:  Number of extruded elements during extrusion. Default is 20. 
- <strong>`element_order`</strong>:  Finite element order. Default is 1 
