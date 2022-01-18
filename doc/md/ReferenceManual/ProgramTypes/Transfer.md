@page transfer_ref Transfer


<strong>Transfer JSON Template</strong>
    
    {
        "Program Type": "Transfer",
        "Mesh File Options": {
            "Source Mesh File": "source"
            "Target Mesh File": "target"
            "Output Mesh File": "output"
        },
        "Transfer Options": {
            ...
        }
    }

\subsection transfer_params Transfer Parameters

- <strong>`Method`</strong>:  Sets the transfer method. Options are `"Consistent Interpolation"`, `"Conservative Surface Transfer"`, and `"Conservative Volume Transfer"`   
- <strong>`Check Transfer Quality`</strong>:  If `true`, performs a quality check.   
- <strong>`Array Names`</strong>:  Provides a list of arrays.   

There are three transfer methods supported:
- `"Consistent Interpolation"`
- `"Conservative Surface Transfer"` (requires <em>IMPACT</em>)
- `"Conservative Volume Transfer"` (requires <em>supermesh</em>)
