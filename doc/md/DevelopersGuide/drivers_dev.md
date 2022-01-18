@page drivers_dev

The driver classes (which inherit from `NEM::MSH::NemDriver`) are the highest level of abstraction. Each driver represents some workflow that can be executed using the virtual `NemDriver::execute()` method. Users specify parameters via getters and setters. Note that the drivers also can be run from the command line using the `nemosysRun` executable and a <em>JSON</em> file, or in Python via `pyNemosys`.

@section implementation Implementing a Driver

Suppose we want to implement some driver `FooDriver`. Here's what the header of `FooDriver` might look like:

~~~{.c}
    class NEMOSYS_EXPORT FooDriver : public NemDriver {
        public:
            using Files = DriverInOutFiles;

            struct NEMOSYS_EXPORT Opts {
                public:
                    explicit Opts(std::vector<int> member1);
                    std::vector<int> member1;
                    double member2{0.9};
                    jsoncons::optional<int> member3{};
                    JSONCONS_TYPE_TRAITS_FRIEND
                private:
                    Opts() = default;
            }

            FooDriver(Files files, Opts opts);

            const Files &getFiles() const;
            void setFiles(Files files);
            const Opts &getOpts() const;
            void setOpts(Opts opts);
            void execute() const override;

            JSONCONS_TYPE_TRAITS_FRIEND

        private:
            FooDriver();
            static constexpr const char *programType = "Foo";
            jsoncons::string_view getProgramType() const override;
            
            Files files_;
            Opts opts_;
    }
~~~

Ideally, drivers only do two things: accept parameters for some workflow, and execute it. To that end, note that `execute` is `const`, because there should be no reason to alter the state of the driver when running. The following sections describe some of the members in more detail.

@subsection privateconstructors Private Constructors

The way that <em>NEMoSys</em> implements <em>JSON</em> deserialization requires that classes have default constructors. However, if the class would not otherwise have a default constructor (because some members need to be set by the user in order for the object to be meaningful), make the default constructor private and use the `JSONCONS_TYPE_TRAITS_FRIEND` macro so that the <em>JSON</em> deserialization can use the default constructor. If the default constructor requires an explicit definition, set the required values to any value--the <em>JSON</em> interface will ensure at run-time that the user provides a value. See @ref deserialization for details on <em>JSON</em> deserialization.


@subsection filesopts Files, files_, Opts, and opts_

Drivers act on files and a set of parameters. To group these, most drivers have nested classes (with standardized names `Files` and `Opts`) and one member of each of these types (note the trailing underscore for private members, following the Google C++ style guide). The `struct` keyword rather than `class` is used to suggest that the purpose of these nested classes is just to hold data. Note that many drivers have just a single output mesh file, or have an input and an output. Use an alias declaration, if necessary (with the `using` keyword) to enforce consistent naming.

Drivers have both required and optional parameters. To enforce this, make sure any publicly available constructors have some way to set each required parameter, so that all `Opts` structs are always in a valid state. Member variables with default values should be initialized where they are declared, rather than through the constructor, both so that it is obvious from the header alone what the default value is, and so that the default constructor can be defined using `= default`. Note that in some cases, options can be turned on or off, but turning them on requires additional parameters. To represent these types of optional values, use the `jsoncons::optional<T>` template, which represents an object of type `T` or the absence of any object at all. For more info about `jsoncons::optional`, note that is has the same interface as `std::optional`.

@subsection deserialization JSON Deserialization for Drivers

Deserialization of a type from <em>JSON</em> relies on macros that generate code defining a specialization of the `jsoncons::json_type_traits<Json, T>` template. Please read [this page](https://github.com/danielaparker/jsoncons/blob/v0.159.0/doc/ref/json_type_traits.md) from the <em>jsoncons</em> documentation. Note that in addition to the macros documented [here]({https://github.com/danielaparker/jsoncons/blob/v0.159.0/doc/ref/json_type_traits/convenience-macros.md), <em>NEMoSys</em> provides additional macros to support virtual inheritance in `NemJsonMacros.H`. Specializations for many STL types (including more containers, `std::pair`, `std::optional`, `std::shared_ptr`, `std::unique_ptr`, and `jsoncons::optional`) are [already provided](https://github.com/danielaparker/jsoncons/blob/v0.159.0/doc/ref/json_type_traits/built-in-specializations.md) based on specializations of `jsoncons::json_type_traits<Json,T>` for the templated type. The following is an example of using these macros for `FooDriver`:

~~~{.c}
    #ifndef NEMOSYS_FOOJSON_H_
    #define NEMOSYS_FOOJSON_H_

    #include "Drivers/FooDriver.H"

    #include "Drivers/NemJsonMacros.H"

    NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
        NEM::DRV::FooDriver, NEM::DRV::NemDriver, 2,
        (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
        (getOpts, setOpts, "Foo Driver Options"),
        (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
            [](const jsoncons::string_view &x) {
                return x == NEM::DRV::FooDriver::programType;
            }))

    JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::FooDriver::Opts, 1,
        (member1, "First Option"),
        (member2, "Second Option"),
        (member3, "Third Option"))

    #endif // NEMOSYS_FOOJSON_H_
~~~

This file might be named `FooJson.H` and live under `include/Drivers` (possibly within a subdirectory, mimicking where `FooJson.C` is located). Note that `FooJson.H` should not be installed! Finally, make sure to add `FooDriver` to `DriverJsonTypeTraits.H`:

~~~{.c}
     NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_BASE(
        NEM::DRV::NemDriver,
        (... // other drivers
            NEM::DRV::FooDriver),
        1,
        (getProgramType, , NEM::DRV::JSON::programType, JSONCONS_RDONLY,
            NEM_JSON_CHECK_KEY_ONLY))
    ...
    #include "Drivers/FooJson.H"
~~~

Only `DriverJsonTypeTraits.H` should be included in any other source file (and not in any public header). Note that the <em>NEMoSys</em> drivers <em>JSON</em> serializations have some standardized names and members. The classes that directly inherit from `NEM::DRV::NemDriver` should have a read-only string-valued key called "Program Type" in the <em>JSON</em> serialization. Note the pure virtual function `getProgramType` in `NEM::DRV::NemDriver` to enforce this. These classes should (1) have a `static constexpr const char *` member called `programType`, (2) implement `getProgramType` by returning `programType`, and (3) check that the "Program Type" in the <em>JSON</em> serialization matches `programType`. To make sure that the serialized name is the same across all drivers, note that `DriverJsonTypeTraits.H` defines the `static constexpr NEM::DRV::JSON::programType` so we can avoid typos in repeating the string "Program Type". Use the `NEM::DRV::JSON` namespace for repeated serialized names!

@subsection jsonconsmacros JSONCONS_ macros

The macros that have the prefix `JSONCONS` come directly from the <em>jsoncons</em>library and are well documented [here](https://github.com/danielaparker/jsoncons/blob/v0.159.0/doc/ref/json_type_traits/convenience-macros.md). The most important parts are the table describing the naming conventions and the table describing the parameters of the macros.

@subsection nemjsonmacros NEM_JSON_ macros

The following additional macros are defined in `NemJsonMacros.H`:

~~~{.c}
    NEM_JSON_N_MEMBER_NAME_TRAITS_BASE(class_name,
        (child_class, child_class, ...),
        num_mandatory,
        (member0,serialized_name0[,mode0,match0,into0,from0]),
        (member1,serialized_name1[,mode1,match1,into1,from1])...)
        
    NEM_JSON_N_MEMBER_NAME_TRAITS_INTERMEDIATE(class_name,
        (child_class, child_class, ...),
        parent_class,
        num_mandatory,
        (member0,serialized_name0[,mode0,match0,into0,from0]),
        (member1,serialized_name1[,mode1,match1,into1,from1])...)
        
    NEM_JSON_N_MEMBER_NAME_TRAITS_FINAL(class_name,
        parent_class,
        num_mandatory,
        (member0,serialized_name0[,mode0,match0,into0,from0]),
        (member1,serialized_name1[,mode1,match1,into1,from1])...)
        
    NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_BASE(class_name,
        (child_class, child_class, ...),
        num_mandatory,
        (getter0,setter0,serialized_name0[,mode0,match0,into0,from0]),
        (getter1,setter1,serialized_name1[,mode1,match1,into1,from1])...)
        
    NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_INTERMEDIATE(class_name,
        (child_class, child_class, ...),
        parent_class,
        num_mandatory,
        (getter0,setter0,serialized_name0[,mode0,match0,into0,from0]),
        (getter1,setter1,serialized_name1[,mode1,match1,into1,from1])...)
        
    NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(class_name,
        parent_class,
        num_mandatory,
        (getter0,setter0,serialized_name0[,mode0,match0,into0,from0]),
        (getter1,setter1,serialized_name1[,mode1,match1,into1,from1])...)
~~~

These macros provide more support for polymorphic classes. They operate the same way (and follow the same naming convention) as the `JSONCONS_` macros with the following additions:

@subsubsection ptrs Pointers

The macros above actually provide partial specializations for `jsoncons::json_type_traits<Json, class_name *>` as opposed to `jsoncons::json_type_traits<Json, class_name>`. Thus, `jsoncons::json_type_traits<Json, class_name *>::as` results in a dynamically allocated object for which the caller takes ownership. Note that `DriverJsonTypeTraits.H` already defines the `jsoncons::json_type_traits<std::unique_ptr<NEM::DRV::NemDriver>>` and `jsoncons::json_type_traits<std::shared_ptr<NEM::DRV::NemDriver>>` specializations.  

@subsubsection polymorphism Polymorphism

`BASE`, `INTERMEDIATE`, and `FINAL` represent where a class is in a class hierarchy:
- `BASE` refers to classes with no base classes (more precisely, no base classes with members that need to be serialized/deserialized) but with some derived classes.
- `INTERMEDIATE` refers to classes with both base and derived classes.
- `FINAL` refers to classes that have no derived classes but some base class (that is, effectively `final`).

Deserialization for types that use these macros will first use `jsoncons::json_type_traits<>::is` to determine the most derived class (the macros use the `num_mandatory`, `serialized_name`, and `match` parameters to construct the `is` methods), then deserialize the members defined in the macro of the most-derived class, then proceed back up the class hierarchy to deserialize remaining members. In particular, this means that, when using these macros, members that are already listed in a base class's macros do not need to be repeated in the derived class (unless they need to be treated differently--see the paragraph on `modeN` extensions).

@subsubsection typevalcheck Checking Type and Value
The `JSONCONS_` macros can check for the presence of mandatory serialized names of a member based on the `num_mandatory` parameter alone. Keep in mind that mandatory names must come first, followed by optional names. They can also check that the corresponding value is valid using the `matchN` parameter. However, there are two pitfalls. First, even if the key would not be checked based on the ordering of the data member and the `num_mandatory` parameter, it <em>is</em> checked if `matchN` is present, so there is no way to define `matchN` for an optional data member. Second, suppose all we want to verify is that some data member is an array of strings as opposed to a single string. Using `matchN` would require actually creating the data object, whereas a call to `jsoncons::json_type_traits<std::array<std::string>>::is` would be sufficient and cheaper.
To resolve these issues, the `NEM_JSON_` macros have slightly different behavior. First, the value of a data member is checked only if the serialized name is present. Second, if the `modeN` parameter is present and the `matchN` parameter is not, then the `jsoncons::json_type_traits<Json, class_name>::is` will call `jsoncons::json_type_traits<Json, member_type>::is` to check the value. If you need to specify `modeN`, `matchN`, `intoN`, or `fromN` but do not want to check the value of the member in `jsoncons::json_type_traits<Json, class_name *>::is}`, use the 
`NEM_JSON_CHECK_KEY_ONLY` macro as the `matchN` parameter.

@subsection othermacros Other NEM_JSON_ macros

There are a few more, less commonly used, macros defined in `NemJsonMacros.H`:

~~~{.c}
NEM_JSON_WRAP_SMART_PTR(class_name, pointer_template)
NEM_JSON_SMART_PTR_VAL(class_name, pointer_template)
NEM_JSON_N_MEMBER_NAME_TRAITS_VAL(class_name,num_mandatory,
    (member0,serialized_name0[,mode0,match0,into0,from0]),
    (member1,serialized_name1[,mode1,match1,into1,from1])...)
~~~


`NEM_JSON_WRAP_SMART_PTR` defines `jsoncons::json_type_traits<Json, pointer_template<class_name>>` using
the template specialization `jsoncons::json_type_traits<Json, class_name *>` (which are what the main `NEM_JSON_` macros above define).  
`NEM_JSON_SMART_PTR_VAL` defines `jsoncons::json_type_traits<Json, pointer_template<class_name>>` using the template specialization
`jsoncons::json_type_traits<Json, class_name>`. Note that the specializations for *non-polymorphic* types are [already provided](https://github.com/danielaparker/jsoncons/blob/v0.159.0/doc/ref/json_type_traits/built-in-specializations.md)  
`NEM_JSON_N_MEMBER_NAME_TRAITS_VAL` works exactly as `JSONCONS_N_MEMBER_NAME_TRAITS`, except it additionally supports the `modeN` extensions described above. If you are invoking this template, however, that probably means the class interface needs to be redesigned!

@section pythonbindings Python Bindings for Drivers

The Python bindings for <em>NEMoSys</em> are defined in `python/pyNemosys.i` and are constructed using <em>SWIG</em> (specifically, <em>SWIG</em> 3.0.12). The `python/CMakeLists.txt` configures a `setup.py` using `setup.py.in`. The `python/CMakeLists.txt` then builds `pyNemosys` by invoking `python setup.py build_ext build`. For details on the `build_ext` command, see the source [here](https://github.com/pypa/setuptools/blob/main/setuptools/command/build_ext.py) which relies on [this](https://github.com/pypa/setuptools/blob/main/setuptools/_distutils/command/build_ext.py). <em>SWIG</em> generates two files: `pyNemosys_wrap.cpp`, the source for a Python extension module named `_pyNemosys`, and `pyNemosys.py`, a Python module that defines proxy Python classes for <em>NEMoSys</em> classes using methods defined in `_pyNemosys`. For details on <em>SWIG</em>, see the [documentation](http://www.swig.org/Doc3.0/SWIGDocumentation.html). Note that there are a few features that are not covered in the documentation; in particular, [this file](https://github.com/swig/swig/blob/v3.0.12/Lib/python/pyuserdir.swg) contains a
few undocumented <em>SWIG</em> directives. The following additions to `pyNemosys.i` show how `FooDriver` might be wrapped:

~~~{.c}
#include <Drivers/FooDriver.H>
%}
%include <Drivers/FooDriver.H>
EXPOSE_FLATNESTED_PY_CLASSES(FooDriver, Opts)
%pythoncode {
FooDriver.Files = DriverInOutFiles
};
CUSTOM_PROPERTIES_INNER_NS(NEM::DRV, FooDriver, Opts, member1)
~~~

@subsection include include

Note the two `include`s above serve two different purposes. The first `include` instructs <em>SWIG</em> to make sure that `pyNemosys_wrap.cpp` includes `Drivers/FooDriver.H` (the purpose of the `%{ %}` block is to include code verbatim in `pyNemosys_wrap.cpp`). The second `include` is a preprocessor directive for the current file. Thus, after <em>SWIG</em> preprocesses `pyNemosys.i`, it will wrap any  classes, functions, and global variables in the preprocessed file, including `FooDriver`.

@subsection namesandnests Namespaces and Nested Classes

By default, <em>SWIG</em> will ignore namespaces when choosing the Python name for classes, so `NEM::DRV::FooDriver` becomes `FooDriver` (note, however, that the Python classes all live inside the `pyNemosys` module). Because the `flatnested` feature is turned on, <em>SWIG</em> will wrap nested classes, but not as attributes of the enclosing class. By default, nested classes will be
attributes of the `pyNemosys` module. However, this causes collision of type names. `pyNemosys.i` renames these types by prepending an
underscore, the enclosing type name, and another underscore. Note the leading underscore, following the Python convention, to suggest internal use.
Thus, `NEM::DRV::FooDriver::Opts` will be wrapped as `_FooDriver_Opts`. Finally, using the `EXPOSE_FLATNESTED_PY_CLASSES` macro, we add the wrapped class as an attribute of the outer class with the expected name so that `_FooDriver_Opts` can be accessed as `FooDriver.Opts` in Python.

@subsection typealiases Type Aliases

<em>SWIG</em> correctly parses type aliases when generating wrappers, but does not provide the analogous Python attributes. Thus we must manually add attributes using the `%pythoncode` directive, which copies code into `pyNemosys.py`.

@subsection stlcontainertypes STL Container Types

See [this section](http://www.swig.org/Doc3.0/SWIGDocumentation.html#Library_stl_cpp_library) of the <em>SWIG</em> documentation on how to wrap STL container types. Note that for types without public default constructors, you must manually ignore constructors that take a size and resize methods.

@subsection containersasdata Containers as Data Members

By default, in order to enforce type checking, <em>SWIG</em> wraps public data members as Python [property](https://docs.python.org/3/library/functions.html#property) attributes. They can be set like:

~~~{.py}
obj = pyNemosys.FooDriver.Opts([0, 1])
obj.member2 = 0
obj.member1 = pyNemosys.VecInt([1, 2])
~~~

Note how the constructor converts automatically between the Python list and `std::vector<int>`. The default behavior for member setters does not allow this, instead requiring the user to explicitly create a `std::vector` before assigning to the data member. There is a workaround, however: use the `CUSTOM_PROPERTIES_` macros and list all container members. Then, assignments like the following will work:

~~~{.py}
obj = pyNemosys.FooDriver.Opts([0, 1])
obj.member1 = [1, 2]
~~~

How does this work? The getter and setter methods that <em>SWIG</em> generates in `pyNemosys_wrap.cpp` for a member of type `T` return and take `T *`.
This is appropriate for the getters; however, for setters, <em>SWIG</em> won’t do automatic conversions from Python types (except the corresponding Python proxy class), just like temporary values can’t bind to lvalue references in C++. Ṫhe workaround is to extend the class with a getter that returns `T &` and a setter that takes `T`, and monkey-patching the data member as a Python `property`.

@subsection optionaljsoncons

Just like the STL containers, we need to make <em>SWIG</em> aware of every `jsoncons::optional<>` type that we plan on using. Use the `OPTIONAL_TEMPLATE_IMMUTABLE` macro for types that are immutable in the Python interface, and `OPTIONAL_TEMPLATE` otherwise (`OPTIONAL_TEMPLATE_IMMUTABLE` turns on `naturalvar` for the type, affecting members of that type). Note that the Python `None` object corresponds to not having a value (the default-constructed `jsoncons::optional`).
