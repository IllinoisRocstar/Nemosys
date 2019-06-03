#include <MeshQualityDriver.H>
#include <MeshQuality.H>
#include <memory>

#ifdef HAVE_CFMSH 
#include <cfmeshQualityParams.H>
#endif

MeshQualityDriver::MeshQualityDriver(std::string _mesh, std::string ofname)
{
    // default constructor does standard check mesh process
    // no improvmenet should be expected
    mesh = meshBase::Create(_mesh);
    mesh->checkMesh(ofname);
    std::cout << "MeshQualityDriver created" << std::endl;
}

MeshQualityDriver::~MeshQualityDriver()
{
    if (mesh)
    {
      delete mesh;
      mesh = 0;
    }
    std::cout << "MeshQualityDriver destroyed" << std::endl;
}

MeshQualityDriver* MeshQualityDriver::readJSON(json inputjson)
{
    MeshQualityDriver* qualdrvobj;
    std::string _mesh;
    std::string ofname;
    std::string engine;
    _mesh = inputjson["Input Mesh File"].as<std::string>();
    ofname = inputjson["Output File"].as<std::string>();
    engine = inputjson.get_with_default("Mesh Quality Engine","default");
    
    if (!inputjson.has_key("Schedule") || (engine == "default"))
    {
        // perform a simple check mesh
        qualdrvobj = new MeshQualityDriver(_mesh, ofname); 
        return qualdrvobj;
    }

    if (!engine.compare("cfmesh"))
    {
        #ifndef HAVE_CFMSH
          std::cerr << "Nemosys must be recompiled with cfMesh support" << std::endl;
          exit(1);
        #endif
    }

#ifdef HAVE_CFMSH 
    cfmshQualityParams* params = new cfmshQualityParams();
    MeshQuality* mq = new MeshQuality(params);
  
    // carry out schedules
    for (auto jsch : inputjson["Schedule"].array_range())
    {
        std::string qImpMtd = jsch["Method"].as<std::string>(); 
        if (!qImpMtd.compare("meshOptimizer"))
        {
            params->nIterations = jsch["Params"]["NIterations"].as<int>();
            params->nLoops = jsch["Params"]["NLoops"].as<int>();
            params->qualThrsh = jsch["Params"]["QualityThreshold"].as<double>();
            params->nSrfItr = jsch["Params"]["NSurfaceIterations"].as<int>();
            if (jsch["Params"].has_key("ConstrainedCellSet"))
            {
                params->_withConstraint = true;
                params->consCellSet = 
                    jsch["Params"]["ConstrainedCellSet"].as<std::string>();
            }

            mq->cfmOptimize();
        }
    }

    return NULL;
#endif
    
}
