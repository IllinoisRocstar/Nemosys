@page geoMeshBase_dev geoMeshBase

`geoMeshBase` is an abstract base class representing a mesh and an optional geometry.

@section memorymanagement Memory Management

`geoMeshBase` inherits from `vtkDataObject`, which takes care of reference counting if used properly. Users of a `fooGeoMesh` class deriving from `geoMeshBase` should create a new `fooGeoMesh` instance in one of the following wats:

~~~{.c}
    void doNothing(fooGeoMesh *fgm) {}

    vtkNew<fooGeoMesh> fgm1;
    doNothing(fgm1);

    // vtkSmartPointer<fooGeoMesh> fgm2 = fooGeoMesh::New() is a mistake
    // because the New() on the right and the assignment both
    // increment the reference count
    vtkSmartPointer<fooGeoMesh> fgm2 = vtkSmartPointer<fooGeoMesh>::New();
    doNothing(fgm2);

    // Also see vtkSmartPointer<>::TakeReference
    vtkSmartPointer<fooGeoMesh> fgm3 =
        vtkSmartPointer<fooGeoMesh>::Take(fooGeoMesh::New());
    doNothing(fgm3);

    fooGeoMesh *fgm5 = fooGeoMesh::New();
    doNothing(fgm5);
    fgm5->Delete();
~~~
@warning Note that using a raw pointer requires a call to `Delete()` and not the C++ keyword `delete` because other objects may still hold references; `Delete()` calls `delete` if the reference count reaches 0. `vtkSmartPointer<>` or `vtkNew<>` are preferred to avoid manual `Delete()` calls.

Note that none of these are on the stack! <em>VTK</em> hides constructors and disables moves and copies.

@section geomesh GeoMesh

`geoMeshBase` objects store data inside the `GeoMesh` type `_geoMesh` data member. Let's take a look at the `GeoMesh` struct (note that it is a protected inner class of `geoMeshBase`).
~~~{.c}
    struct GeoMesh {
        vtkSmartPointer<vtkUnstructuredGrid> mesh;
        std::string geo;
        std::string link;
        SideSet sideSet;
    };
~~~
The four members of `GeoMesh` are:
    - `mesh` contains the mesh. It should never be `nullptr`. If a `geoMeshBase` object needs to be reset for some reason, set this field with `vtkSmartPointer<vtkUnstructuredGrid>::New()`.
    - `geo`, if not empty, is the name of a gmsh model representing the geometry. Remember to call `gmsh::model::setCurrent` before other <em>Gmsh</em> API calls.
    - `link`, if it and `geo` are not empty, is the name of a `vtkDataArray` in `mesh->GetCellData()` that denotes, for every cell in `mesh`, which physical group of `geo` it belongs to.
    - `sideSet` represents the edges (if the mesh is 2D) or faces (if the mesh is 3D) on the boundaries of physical groups. Therefore, every cell in `sideSet.sides` ("side cell") is the edge/face of some cell in `mesh` ("original cell"). Note that `sideSet.sides->GetPoints()` and `mesh->GetPoints()`, if they both exist, should be pointers to the same object, so that point indices remain consistent. The `sideSet.sides` has one required cell data array (`"GeoEnt"`), two arrays that are required if the mesh type represents side sets as sides of a cell (as opposed to <em>Gmsh</em> for example, which represents the mesh elements of Physical Surfaces independently form the elements on any Physical Volume) (`"OrigCellIDs"` and `"CellFaceIds"`), and an optional cell data array (`"TwinIds"`):
        - `"GeoEnt"`, aka `GeoMeshBase::SIDE_SET_GEO_ENT_NAME`, is a `vtkIntArray` which, for each side cell, contains the physical group that it belongs to in `geo`.
        - `"OrigCellIDs"` is a `vtkIdTypeArray` which, for each side cell, contains the index in `mesh` of the original cell. 
        - `"CellFaceIds"` is a `vtkIntArray` which, for each side cell, contains the edge/face of the original cell that has the same set of points (ordering not guaranteed, although normals should have the same direction). Note that this follows <em>VTK</em> indexing, so that an edge/face can be recovered from `vtkCell::GetEdge` and `vtkCell::GetFace`, respectively.
        - `"TwinIds"` is a `vtkIdTypeArray` which, for each side cell, contains the index in `sideSet.sides` of the twin side across a material interface (same points, opposite normal) or `-1` if not such twin exists.

Use the constructor/getters/setters for these arrays so that the array names are set correctly. Subclasses of `geoMeshBase` can use other cell/point/field data arrays if necessary.

@section implementation Implementing a geoMeshBase class

Suppose for some mesh type `Foo`, we want to implement a `fooGeoMesh` class that holds a `Foo` mesh. Here's what the header of the class might look like. Be sure to check ~FILL IN SECTION~ for necessary changes outside of the `fooGeoMesh` class.
~~~{.c}
    class fooGeoMesh : public geoMeshBase {
        public:
            vtkTypeMacro(fooGeoMesh, geoMeshBase)
            static fooGeoMesh *New();
            static fooGeoMesh *Read(const std::string &fileName);
            void write(const std::string &fileName) override;
            void report(std::stream &out) const override;
            const Foo &getFooMesh() const;
            void setFooMesh(Foo fooMesh);
        protected:
            fooGeoMesh();
            explicit fooGeoMesh(Foo fooMesh);
        private:
            static GeoMesh Foo2GM(const Foo &fooMesh);
            static Foo GM2Foo(const GeoMesh &geoMesh);
            void resetNative() override;
            Foo fooMesh_{};
    }
~~~
Note that `fooGeoMesh` has a `Foo` type data member. `geoMeshBase` is not designed to be aware of external changes to the native mesh type, so `geoMeshBase` subclasses should own (whether as a data member or using an owning pointer) the mesh of the native type.

The next sections discuss the details of the methods in the header (and other necessary changes!) and tips on implementing them, given in order of suggested implementation. Note that while `Foo2GM` and `GM2Foo` are technically optional, implementing these methods is highly recommended to make implementation of other methods easier.

@subsection vtktypemacro vtkTypeMacro

The `vtkTypeMacro` declares and defines methods that <em>VTK</em> uses to do run-time type checking and cloning. The macro generates both declarations and definitions.

@subsection new New

The static `New` method is used to create new instances and is also how `vtkSmartPointer<>` and `vtkNew<>` create new instances.
Sample implementation:
~~~{.c}
    vtkStandardNewMacro(fooGeoMesh)
~~~
@subsection foo2gm Foo2GM

`Foo2GM` should construct a `GeoMesh` struct that contains the information in `fooMesh_`. Its purpose is to simplify the implementation of the `fooGeoMesh` (`FooMesh`) constructor (because it should call the `geoMeshBase` (`GeoMesh`) constructor in the member initializer list, which can't declare any variables). Recall that the `GeoMesh` struct has four members (see @ref geomesh ). If `fooMesh_` is empty of invalid, set the `mesh` member to be empty (using `vtkSmartPointer<vtkUnstructuredGrid>::New()`) as opposed to `nullptr`. Only set information if it is readily available in `fooMesh_`! If a `Foo` has no concept of a geometry, leave `geo` and `link` as empty strings. If finding the boundary sides is expensive, use the default-constructed `sideSet`.
Sample implementation:
~~~{.c}
    static geoMeshBase::GeoMesh fooGeoMesh::Foo2GM(const Foo &fooMesh) {
        auto vug = vtkSmartPointer<vtkUnstructuredGrid>::New();
        std::string geo{};
        std::string link{};
        SideSet sideSet{};
        vug->SetPoints(...);
        vug->InsertNextCell(...);
        if (fooMesh.hasGeometry()) {
            GmshInterface::Initialize();
            std::string geo{"geoMesh_" + nemAux::getRandomString(6)};
            gmsh::model::add(geo);
            gmsh::model::setCurrent(geo);
            std::string link = GEO_ENT_DEFAULT_NAME;
            gmsh::model::mesh::addPhysicalGroup(...);
            auto linkArr = vtkSmartPointer<vtkIntArray>::New();
            linkArr->SetName(link.c_str());
            linkArr->InsertNextValue(...);
            vug->GetCellData()->AddArray(linkArr);
            auto sideSetPD = vtkSmartPointer<vtkPolyData>::New();
            sideSetPD->SetPoints(vug->GetPoints());
            sideSetPD->InsertNextCell(...);
            vtkNew<vtkIntArray> sideSetEntArr;
            sideSet = {sideSetPD, sideSetEntArr};
        }
        return {vug, geo, link, sideSet};
    }
~~~

@subsection gm2foo GM2Foo

`GM2Foo` should construct a `Foo` object from a `GeoMesh` struct. Its purpose is to simplify the implementation of `resetNative`. In particular, using a separate static method reduces the chances of coding mistakes by preventing access to the old `fooMesh_` that should be overwritten in `resetNative`. Note that you should alter the signature as necessary in the case that `Foo` meshes are better represented by more than one object.

@subsection constructors Constructors

The default constructor is necessary for `New`. It should initialize `fooMesh_` and `GeoMesh` to be empty meshes (if using owning pointers to represent the native mesh type, valid objects should be preferred over `nullptr` to represent empty meshes). The easiest way to do this is to delegate to the `fooGeoMesh(Foo)` constructor, which itself should use a member initializer list to call the `geoMeshBase(GeoMesh)` constructor, for example:
~~~{.c}
    fooGeoMesh::fooGeoMesh() : fooGeoMesh(Foo{}) {}
    fooGeoMesh::fooGeoMesh(Foo fooMesh)
    : geoMeshBase(Foo2GM(fooMesh)), fooMesh_(std::move(fooMesh)) {}
~~~

@subsection resetnative resetNative

This method should overwrite all `Foo`-specific data members of `fooGeoMesh` by using the `GM2Foo` method on `getGeoMesh`. Note that this may require calling `geoMeshBase::findSide2OrigCell` if `Foo` requires the original cell/face for elements in the `geoMeshBase::sideSet`. The `resetNative` method is used by `takeGeoMesh`, the canonical way of converting between `geoMeshBase` types, so `Foo`-specific data members should be treated as stale.
Sample implementation:
~~~{.c}
    void fooGeoMesh::resetNative() {
        fooMesh_ = GM2Foo(getGeoMesh());
    }
~~~

@subsection getfoomesh getFooMesh

This should return a view/read-only version of `fooMesh_`. Ideally the return type should be `const Foo &`. If this is not possible, return a copy of `fooMesh_`. Again, this is because `geoMeshBase` is not aware of external changes to `fooMesh_`.
Sample implementation:
~~~{.c}
    const FooMesh &fooGeoMesh::getFooMesh() {
        return fooMesh_;
    }
~~~

@subsection setfoomesh setFooMesh

This should reset the object to a state similar to constructing a new object using the `fooGeoMesh(Foo)` constructor.
Sample implementation:
~~~{.c}
    void fooGeoMesh::setFooMesh(Foo fooMesh) {
        this->setGeoMesh(Foo2GM(fooMesh));
        this->fooMesh_ = std::move(fooMesh);
    }
~~~

@subsection read Read

The static `Read` method should take a file path, read the `Foo` mesh from the file, and then construct and return a new `fooGeoMesh`.
Sample implementation:
~~~{.c}
    fooGeoMesh *fooGeoMesh::Read(const std::string &fileName) {
        Foo fooMesh{};
        fooMesh.readFromFile(fileName);
        return new fooGeoMesh(std::move(fooMesh));
    }
~~~

@subsection write write

Write the `fooMesh_` to a file. `fooGeoMesh` should support the file types that `Foo` supports. Avoid accessing the `GeoMesh`; prefer accessing `Foo`-specific data only.

@subsection report report

Write a short summary of the mesh to a stream. Feel free to make use of `geoMeshBase::report`.
Sample implementation:
~~~{.c}
    void fooGeoMesh::report(std::ostream &out) const {
        geoMeshBase::report(out);
        out << fooMesh_.summary() << '\n';
    }
~~~

@subsection others Other changes

Be sure to:
- Add an enumerator to `MeshType` that corresponds to `fooGeoMesh`.

~~~{.c}
    enum class MeshType {
        ...,
        FOO_GEO_MESH
    }
~~~

- Add a new case in `MeshTypeFromFilename` for file extensions that are associated with `Foo` objects.

~~~{.c}
    MeshType MeshTypeFromFilename(const std::string &fileName) {
    std::string fileExt = nemAux::find_ext(fileName);
    ...
    } else if (fileExt == ".foo") {
        return MeshType::FOO_GEO_MESH;
    } else {
    ...
}
~~~

- Add cases to `Read` and `New` for the new `MeshType` enumerator.

~~~{.c}
    geoMeshBase *Read(const std::string &fileName, MeshType meshType) {
        switch (meshType) {
            ...
            case MeshType::FOO_GEO_MESH: return fooGeoMesh::Read(fileName);
        }
    }

    geoMeshBase *New(MeshType meshType) {
        switch (meshType) {
            ...
            case MeshType::FOO_GEO_MESH: return fooGeoMesh::New();
        }
    }
~~~

- Add a case for `fooGeoMesh` in `srvBase::RequestDataObject`.

~~~{.c}
    int srvBase::RequestDataObject(vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector) {
        ...
        } else if (typeName == "fooGeoMesh") {
            output = MSH::fooGeoMesh::New();
        } else if (typeName == "geoMeshBase") {
            ...
        }
    }
~~~
