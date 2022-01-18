@page json_anatomy Anatomy of the JSON file

@tableofcontents

The standard input file for <em>Promesh</em> is in JSON format. It typically has the following format:

    {
        "Program Type": "<PROGRAM_TYPE>",
        "Mesh File Options": {
            ...
        },
        "<PROGRAM_TYPE> Options": {
            ...
        }
    }

where `"<PROGRAM_TYPE>"` can be substituted for one of the following:
- `"Verification"`
- `"Conversion"`
- `"Input Generation"`
- `"Mesh Generation"`
- `"Mesh Quality"`
- `"NucMesh Generation"`
- `"Pack Mesh Generation"`
- `"Refinement"`
- `"Template Mesh Generation"`
- `"Transfer"`

Two other options, `"Proteus"` and `"Rocstar Communication Generation"` are not detailed in this guide.

Many program types will also require additional top-level keywords; check the reference for each type to verify that all the necessary keywords are present.

The `"Mesh File Options"` section includes specifications for the files required. This may be geometry input file, mesh output files, or mesh input files. The `"<PROGRAM_TYPE> Options"` section includes all the keywords associated with that program type.


<strong>JSON Syntax</strong>

From [here](https://www.json.org/json-en.html):

JSON (JavaScript Object Notation) is a lightweight data-interchange format. It is easy for humans to read and write. It is easy for machines to parse and generate.

JSON is built on two structures:
- A collection of name/value pairs, called an object
- An ordered list of values, called an array

Objects are enclosed in curly braces. Each name is followed by a colon, and each name/value pair is separated by a comma. Arrays are enclosed in square brackets, and values within arrays are separated by commas.

Values can be strings (in double quotes), numbers, `true` or `false` or `null`, objects, or arrays. These structures can be nested   

