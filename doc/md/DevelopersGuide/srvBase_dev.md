@page srvBase_dev srvBase

`srvBase` is an abstract base class that implements complex operations on `geoMeshBase` objects. Note that `srvBase` inherits from `vtkAlgorithm`, so the memory management tips in @ref geomeshbase also apply to `srvBase` objects. Some additional resources for understanding `vtkAlgorithm`:

- This [wiki page](https://vtk.org/Wiki/VTK/Tutorials/New_Pipeline) provides an overview of `vtkAlgorithm`.
- The [VTK Users' Guide](https://www.kitware.com/products/books/VTKUsersGuide.pdf) chapters 15 and 17 offer a more detailed guide to implementing custom `vtkAlgorithm`s. (Caution: this book was written for VTK 5.4!)
- The new [VTK examples site](https://kitware.github.io/vtk-examples/site/) includes an updated overview (see VTKBook > Chapter 4, in particular sections 4.2-4.4) and full \cpp\ examples (see Cxx > Developers).
- This [series of blog posts](https://blog.kitware.com/a-vtk-pipeline-primer-part-1/) provides examples (mostly Python) of a custom `vtkAlgorithm` as well as some explanation of the design of `vtkAlgorithm`.

@section execution Execution

Note that `srvBase` objects can be pipelined, that is, services can take the other `srvBase` objects as inputs. More precisely, output ports and input ports of different `srvBase` objects can be connected, as in:

~~~{.c}
    vtkNew<fooGeoMesh> mesh;
    vtkNew<myFooSrv> srv1;
    vtkNew<myFooSrv> srv2;
    srv1->SetInputDataObject(mesh);
    srv2->SetInputConnection(srv1->GetOutputPort());
    srv2->Update();
~~~
Execution of the entire pipeline only happens when `Update` is called on `srv2`. If the service has not been modified since the last time it was updated, then it does nothing. Otherwise, the service and any inputs to the service that need to be updated (and their inputs and so on) do four things, in the following order:

- `RequestDataObject` is called on each service in the forward direction, meaning `RequestDataObject` is called on a service's inputs before it is called on the service (`srv1` before `srv2`). This method creates empty objects to store the results of a service.
- `RequestInformation` is called on each service in the forward direction. This method passes metadata about the result of a service to subsequent services that might need it.
- `RequestUpdateExtent` is called on each service in the backward direction, meaning `RequestDataObject` is called on a service before it is called on the service's inputs (`srv2` before `srv1`). This method tells its input services what portion of the data it needs, in case the input service can save resources by only partially executing.
- `RequestData` is called on each service in the forward direction. `RequestData` runs the actual algorithm that a service represents by writing data into the object previously created by `RequestDataObject`.

Note that the input and output information and data are all passed between services via `vtkInformation`, which holds key-value pairs, where the keys are instances of `vtkInformationKey`.

@section implementation Implementing a srvBase class

Suppose we want to implement some service `myFooSrv` that takes `fooGeoMesh` objects, runs some algorithm, and results in another `fooGeoMesh` object, with some option `Bar` that controls the algorithm. Here's what the header of `myFooSrv` might look like:

~~~{.c}
    class myFooSrv : public srvBase {
        public:
            vtkTypeMacro(myFooSrv, srvBase)
            static myFooSrv *New();
            vtkSetMacro(Bar, int);
            vtkGetMacro(Bar, int);
        protected:
            myFooSrv();
            int RequestData(vtkInformation *request,
                vtkInformationVector **inputVector,
                vtkInformationVector *outputVector) override;
            int FillInputPortInformation(int port, vtkInformation *info) override;
            int FillOutputPortInformation(int port, vtkInformation *info) override;
        private:
            int Bar{};
    }
~~~
Note the use of Pascal/upper camel case for the private data member `Bar` so that the setters and getters created by `vtkSetMacro` and `vtkGetMacro` are consistently capitalized.

@subsection vtktypemacro vtkTypeMacro

The `vtkTypeMacro` declares and defines methods that <em>VTK</em> uses to do run-time type checking and cloning. The macro generates both declarations and definitions.

@subsection new New

The static `New` method is used to create new instances, and is also how `vtkSmartPointer<>` and `vtkNew<>` create new instances.
Sample implementation:
~~~{.c}
    vtkStandardNewMacro(myFooSrv)
~~~

@subsection setandget vtkSetMacro and vtkGetMacro

Prefer using these macros (and similar ones in `vtkSetGet.h`) for implementing setters and getters of options of services. In particular, if writing custom setter functions (or any functions that alter the execution of the algorithm), remember to call `Modified` before returning.

@subsection portinfo FillInputPortInformation and FillOutputPortInformation

The inputs and outputs of a service (generally `geoMeshBase` objects; simpler inputs can be treated as options) are defined in terms of ports. Input ports should at least set the `INPUT_REQUIRED_DATA_TYPE()` key so that the input data type can be checked. Other keys that might be of interest include `vtkAlgorithm::INPUT_IS_OPTIONAL()` and `vtkAlgorithm::INPUT_IS_REPEATABLE()`. Output ports should at least set the `DATA_TYPE_NAME()` key so that the output data type can be checked and so that the appropriate data type is used in `RequestDataObject`. Note that the value for `DATA_TYPE_NAME()` should be the same for all instances of the service. If you find yourself tempted to set this value depending on options or the input, you most likely need to override `RequestDataObject` and implement the logic there. The return values indicate success (`1`) or failure (`0`).
Sample implementation:
~~~{.c}
    int myFooSrv::FillInputPortInformation(int port, vtkInformation *info) {
        if (port == 0) {
            info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "fooGeoMesh");
            return 1;
        }
        return 0;
    }

    int myFooSrv::FillOutputPortInformation(int port, vtkInformation *info) {
        if (port == 0) {
            info->Set(vtkDataObject::DATA_TYPE_NAME(), "fooGeoMesh");
            return 1;
        }
        return 0;
    }
~~~

@subsection constructor Constructor

The constructor should set the number of input and output ports for the service.
Sample implementation:
~~~{.c}
    myFooSrv::myFooSrv() {
        this->SetNumberOfInputPorts(1);
        this->SetNumberOfOutputPorts(1);
    }
~~~

@subsection requestdata RequestData

This method contains the main logic of the service. Note that the input and output objects are contained in the `inputVector` and `outputVector` arguments. Be sure to return `0` on failure or `1` on success.
Sample implementation:
~~~{.c}
    int myFooSrv::RequestData(vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector) {
        // first index corresponds to input port
        // second index used if input is repeatable
        auto inInfo = inputVector[0]->GetInformationObject(0);
        // index corresponds to output port
        auto outInfo = outputVector->GetInformationObject(0);
        auto input = MSH::fooGeoMesh::SafeDownCast(
            inInfo->Get(vtkDataObject::DATA_OBJECT()));
        auto output = MSH::fooGeoMesh::SafeDownCast(
            outInfo->Get(vtkDataObject::DATA_OBJECT()));
        Foo outputFoo = input->getFooMesh(); // copy
        outputFoo.runAlgorithm(this->Bar);
        output->setFooMesh(std::move(outputFoo));
        return 1;
    }
~~~

@subsection requestdataobject RequestDataObject

This method creates the output data object for the service before anything is executed. `srvBase` contains a default implementation for this method, which uses the `vtkDataObject::DATA_TYPE_NAME()` key. If the value given by this key refers to a concrete type, then a new instance of the concrete type will be created. If the value associated with `DATA_TYPE_NAME()` is `"geoMeshBase"`, then the output object will be created by calling `NewInstance` on the input object of the input port (specifically, the first object on the input port with the same index). If some custom logic is needed for creating the output object, override this method. Ths overridden method should create the output data object and set it as the value for the `vtkDataObject::DATA_OBJECT()` key for each `vtkInformation` in the `outputVector`. Be sure to return `0` on failure or `1` on success Note that each `FillOutputPortInformation` should set the `DATA_TYPE_NAME()` key even in `RequestDataObject` is overridden.

@subsection requestinfo RequestInformation

This method provides metadata about the output that subsequent services might need. The `srvBase` implementation does nothing; override it if needed.

@subsection requestupdateextent RequestUpdateExtent

This method notifies input services about what extent of the input it needs. The `srvBase` implementation sets `vtkStreamingDemandDrivenPipeline::EXACT_EXTENT()` to 1 for all inputs. Override it if needed, particularly for streaming algorithms.