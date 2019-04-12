// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

// the template parameter file for this executable
#include "Parameters.h"
// the parameter file for the current test case
#ifdef _HAVE_PARSER_
 #include "moveItParse.h"
#else
 #include "MyParams.h"
#endif

#include "MAdLib.h"

#include <iostream>
#include <sstream>
#include <math.h>
#include <sys/time.h>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <set>

#ifdef PARALLEL
 #include "mpi.h"
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::stringstream;
using std::ofstream;
using std::vector;
using std::set;

using namespace MAd;

// ----------------------------------------------------------------------
#ifdef PARALLEL
class EmptyExchanger : public MDB_DataExchanger
{
public:
  EmptyExchanger(int _tag): MDB_DataExchanger(_tag) {}
  ~EmptyExchanger() {}

  void * sendData (pEntity pe,    // in
		   int iProcDest, // in
		   int &_size ) {
    _size = 0;
    return NULL;
  }
  void receiveData (pEntity pe,      //in
	            int iProcSender, //in
		    void *buf ) {}
  void deleteExternalData( pEntity pe) const {}
};
#endif


// ----------------------------------------------------------------------
double CPUTime() {
  
#ifdef PARALLEL
  return MPI_Wtime();
#else
  struct timeval tp;
  struct timezone tz;

  gettimeofday(&tp,&tz);

  return ((double) tp.tv_sec +
          (double) ((double) .000001 * (double) tp.tv_usec));
#endif
}

// ----------------------------------------------------------------------
#if defined(__linux) || defined(linux)
void process_mem_usage(double& vm_usage, double& resident_set)
{
  using std::ios_base;
  using std::ifstream;
  using std::string;
  
  vm_usage     = 0.0;
  resident_set = 0.0;
  
  // 'file' stat seems to give the most reliable results
  //
  ifstream stat_stream("/proc/self/stat",ios_base::in);
  
  // dummy vars for leading entries in stat that we don't care about
  //
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;
  
  // the two fields we want
  //
  unsigned long vsize;
  long rss;
  
  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
              >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
              >> utime >> stime >> cutime >> cstime >> priority >> nice
              >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
  
  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  vm_usage     = vsize / 1024.0 / 1024;
  resident_set = rss * page_size_kb / 1024;
}
#endif

// ----------------------------------------------------------------------
void displayMemoryUsage(string step="")
{
  double ram, virt;

#if defined(__linux) || defined(linux)
  process_mem_usage(virt,ram);
  
#ifdef PARALLEL
    double send[2] = {virt,ram};
    double recv[2];
    MPI_Allreduce(send,recv,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    virt  = recv[0];
    ram = recv[1];
#endif

  std::cout << "Memory usage at step \'"<<step<<"\': " 
            << ram  << " Mb (resident), "
            << virt << " Mb (virtual)\n";
#endif
}

// ----------------------------------------------------------------------
void writeSolution(MeshAdapter* ma, int iter, string type) 
{
  stringstream ss;
  string iterStr;  ss << iter;  ss >> iterStr;
  
  if ( !strcmp(type.c_str(),"Pos") || !strcmp(type.c_str(),"MshAndPos") ) {
    string namePos = "result" + iterStr + ".pos";
    ma->writePos(namePos,OD_MEANRATIO);
//     ma->writePos(namePos,OD_SIZEFIELD_MEAN);
  }
  
  if ( !strcmp(type.c_str(),"Msh") || !strcmp(type.c_str(),"MshAndPos") ) {
    string nameMsh = "result" + iterStr + ".msh";
    ma->writeMsh(nameMsh);
  }

// #warning "output distance"
//   string nameDist = "distance" + iterStr;
//   ma->writeDistanceToWalls(nameDist);

// #warning "output curvature"
//   string nameCurv = "curvature" + iterStr;
//   ma->writeVolumicCurvature(nameCurv);
}

// ----------------------------------------------------------------------
void setMobileObject(mobileObject* mob, pMesh mesh, ObjectDef def)
{
  // --- Name ---
  string name = def.name;
  mob->setName(name);

  // --- Geometry ---
  int type = def.geomLevel;
  set<int>::const_iterator gIter = def.geomTags.begin();
  set<int>::const_iterator gLast = def.geomTags.end();
  for (; gIter != gLast; gIter++) {
    for (int iTag = (*gIter); iTag <= (*gIter); iTag++) {
      mob->addGEntity(type,iTag);
    }
  }

  // --- Kinematics ---
  KinematicType kinType = def.kinType;
  if ( kinType == KT_DISPLACEMENT ) {
    mob->setDxKinematics(def.kinExpression);
  }
  else {
    mob->setVKinematics(PARSED, NULL, NULL, def.kinExpression);
  }

  // --- Local size fields ---
  std::list<LocalSFDef>::const_iterator sIter = def.sizes.begin();
  std::list<LocalSFDef>::const_iterator sLast = def.sizes.end();
  for (; sIter != sLast; sIter++) {
    LocalSizeField* locSField = new LocalSizeField(mesh, (*sIter).name, 
                                                   (*sIter).distToFaces);
    bool   iso = (*sIter).isotropic;
    double radius = (*sIter).radius;
    string sizeN = (*sIter).sizeN;
    string sizeT = (*sIter).sizeT;
    bool   limit = (*sIter).limit;
    double limitSizeTg = (*sIter).limitSizeTg;
    double maxCurv = (*sIter).maxCurv;
    if ( iso ) locSField->setIsoSize(radius,sizeN);
    else       locSField->setAnisoSize(radius,sizeN,sizeT);
    if ( limit ) locSField->setCurvatureLimiter(limitSizeTg, maxCurv);
    
    set<int>::const_iterator gIter = def.geomTags.begin();
    set<int>::const_iterator gLast = def.geomTags.end();
    for (; gIter != gLast; gIter++) {
      for (int iTag = (*gIter); iTag <= (*gIter); iTag++) {
        locSField->addGeometricEntity(type,iTag);
      }
    }
    locSField->updateTree();
    mob->addLocalSField(locSField);

// #warning "debug"
//     locSField->printPosAnisotropic(mesh,"locSF");
  }
}

// ----------------------------------------------------------------------
void deleteObjects(mobileObjectSet* objSet)
{
  set<mobileObject*> objs = objSet->getObjects();
  set<mobileObject*>::iterator oIter = objs.begin();
  for (; oIter != objs.end(); oIter++) {
    
    mobileObject * mob = *oIter;

    // clean size fields
    set<LocalSizeField* > localSFs = mob->getSizes();
    set<LocalSizeField* >::iterator sIter = localSFs.begin();
    for (; sIter != localSFs.end(); sIter++) {
      if (*sIter) delete (*sIter);
    }
    
    // delete object
    if (mob) delete mob;
  }
}

// ----------------------------------------------------------------------
int main(int argc, char* argv[]) 
{
  double main_t0 = CPUTime();

  MAdLibInitialize(&argc,&argv);

  displayMemoryUsage("initialized");

#ifdef _HAVE_PARSER_
  parseControlFile(argc,argv);
#endif

  MoveParameters parameters = getParameters();

  // ------------------------------------------------
  // setup the output
  // ------------------------------------------------

  string outputType = parameters.outType;
  int outputFrequency = parameters.outFrequency;
  string outputPrefix = parameters.outPrefix;

  // ------------------------------------------------
  // load the mesh
  // ------------------------------------------------

  cout << "Loading the mesh...\n";
  double cpu_mesh_0 = CPUTime();

  // --- Reading model ---
  string meshFileName = parameters.meshFileName;
  string geoFileName  = parameters.geoFileName;
  pGModel model = 0;
  GM_create(&model,"theModel");
  if ( !geoFileName.empty() ) GM_read(model, geoFileName.c_str());
  else                        GM_readFromMSH(model, meshFileName.c_str());

  // --- Reading mesh ---
  pMesh mesh = M_new(model);
  M_load(mesh,meshFileName.c_str());
//   M_writeMsh (mesh, (outputPrefix + "initMesh.msh").c_str(), 2, NULL);

  double cpu_mesh_tot = CPUTime() - cpu_mesh_0;
  cout << "Loaded the mesh in "<<cpu_mesh_tot<<" seconds\n";
  displayMemoryUsage("Mesh loaded");

  // ------------------------------------------------
  // build the mesh adapter
  // ------------------------------------------------

  cout << "Building the mesh adapter...\n";
  double cpu_ma_0 = CPUTime();
  MeshAdapter* ma = new MeshAdapter(mesh);

  ma->setOutputPrefix(outputPrefix);

  ma->setMaxIterationsNumber(parameters.maxInnerIter);

  double lowerLenghtBound = parameters.lowerLength;
  double upperLenghtBound = parameters.upperLength;
  ma->setEdgeLenSqBounds( lowerLenghtBound * lowerLenghtBound,
                          upperLenghtBound * upperLenghtBound );

  ma->setSwapMinImproveRatio (parameters.swapMinImproveRatio);
  ma->setSliverQuality       (parameters.sliverQuality);
  ma->setSliverPermissionInESplit    (parameters.splitCanMakeSlivers,
                                      parameters.makeSliverInSplitLenSqBound);
  ma->setSliverPermissionInECollapse (parameters.collapseCanMakeSlivers,
                                      parameters.makeSliverInCollapseLenSqBound);
  ma->setCollapseOnBoundary  (parameters.collapseOnBoundary,
                              parameters.clpOnBdryTolerance);
  ma->setSwapOnBoundary      (parameters.swapOnBoundary,
                              parameters.swapOnBdryTolerance);
  ma->setGeoTracking         (parameters.trackGeometry,
                              parameters.snap_cavityIsMesh,
                              parameters.snap_thickness,
                              parameters.snap_chi,
                              parameters.snap_check,
                              parameters.snap_force);
  ma->setSizeFieldSmoothing  (parameters.SFSmoothing,
                              parameters.SFSmoothGrad);
  ma->setInfiniteLength(parameters.infLength);

// #ifdef PARALLEL
//   EmptyExchanger * dExch = new EmptyExchanger(4238458);
//   ma->setDataExchanger( dExch );
// #endif

  ma->setSFUpdateFrequency(0);

  double cpu_ma_tot = CPUTime() - cpu_ma_0;
  cout << "Built the mesh adapter in "<<cpu_ma_tot<<" seconds\n";

  displayMemoryUsage("Mesh adapter built");

// #warning "debug"
//   ma->writePos("testCurvMaxVec.pos",OD_CURVATURE_MAX_VEC);
//   printf("Written!");

  // ------------------------------------------------
  // build the size fields
  // ------------------------------------------------
  cout << "Buidling the size fields...\n";
  double cpu_sf_0 = CPUTime();

  set<pSField> sizeFields;

  std::list<SizeFieldDef>::const_iterator sIter = parameters.sizes.begin();
  std::list<SizeFieldDef>::const_iterator sLast = parameters.sizes.end();
  for (; sIter != sLast; sIter++) {

    SizeFieldBase* sizeField = NULL;
  
    if ( (*sIter).type == ANALYTICALSFIELD )
      {
        OrientationType orient = (*sIter).orientation;
        if ( orient == ORT_ISOTROPIC ) {
          sizeField = new AnalyticalSField((*sIter).isoSize);
        }
        else {
          sizeField = new AnalyticalSField((*sIter).anisoSize,
                                           (*sIter).anisoDir0,
                                           (*sIter).anisoDir1,
                                           (*sIter).anisoDir2);
        }
      }
    else if ( (*sIter).type == DISCRETESFIELD )
      {
        sizeField = new PWLSField(mesh);
        string source = (*sIter).pwlSource;
        if ( !strcmp(source.c_str(),"InitialLength") ) {
          ((PWLSField*) sizeField)->setCurrentSize();
        }
        else if ( !strcmp(source.c_str(),"Curvature") ) {
          bool aniso = (*sIter).curv_aniso;
          double alpha = (*sIter).curv_alpha;
          double hMin = (*sIter).curv_hMin;
          ((PWLSField*) sizeField)->setCurvatureSize(aniso,alpha,hMin);
        }
        else throw;
      }
    else if ( (*sIter).type == BACKGROUNDSFIELD )
      {
        sizeField = new BackgroundSF();
        string mshFile = (*sIter).mshFile;
        ((BackgroundSF*) sizeField)->loadData(mshFile);
      }
    else throw;

    sizeFields.insert(sizeField);
    ma->addSizeField(sizeField);
  }

  double cpu_sf_tot = CPUTime() - cpu_sf_0;
  cout << "Built the size field in "<<cpu_sf_tot<<" seconds\n";

  displayMemoryUsage("Size fields built");

// #warning "debug"
//   ma->writePos("initSF.pos",OD_SIZEFIELD_MEAN);

  ma->printStatistics(cout);


  // ------------------------------------------------
  // setup the statistics monitoring
  // ------------------------------------------------

  string statFileName = outputPrefix + "modifications";
  FILE* statFile = fopen(statFileName.c_str(),"w");
  if (!statFile) {
    cerr << "Could not open the statistics file " << statFileName << endl;
    throw;
  }

  fprintf(statFile,"iter\t# edge splits\t# edge collaspes\t# edge swaps\t");
  fprintf(statFile,"CPU split\tCPU collapse\tCPU swap\tCPU sliver\n");
  fclose(statFile);


  // ------------------------------------------------
  // setup the quality monitoring
  // ------------------------------------------------

  string qualityFileName = outputPrefix + "quality";
  FILE* qualFile = fopen(qualityFileName.c_str(),"w");
  if (!qualFile) {
    cerr << "Could not open the quality file " << qualityFileName << endl;
    throw;
  }

  double meanShape, worstShape;
  ma->getStatistics(&meanShape,&worstShape);
  fprintf(qualFile,"iter\tmean quality\tworst quality\n");
  fprintf(qualFile,"%d\t%g\t%g\n",0,meanShape,worstShape);
  fclose(qualFile);

  // ------------------------------------------------
  // setup the cpu time monitoring
  // ------------------------------------------------

  string cpuFileName = outputPrefix + "cpu";
  FILE* cpuFile = fopen(cpuFileName.c_str(),"w");
  if (!cpuFile) {
    cerr << "Could not open the cpu report file " << cpuFileName << endl;
    throw;
  }

  fprintf(cpuFile,"iter\tnodes motion\tadaptation\toutputs  \taccumulated nm\taccum adapt\taccum out\ttotal\n");
  fclose(cpuFile);

  string produceCpuInfo = "cat /proc/cpuinfo > \"" + outputPrefix + "cpuinfo\"";
  system(produceCpuInfo.c_str());

  // ------------------------------------------------
  // setup the mesh size monitoring
  // ------------------------------------------------

  string MSFileName = outputPrefix + "meshSize";
  FILE* MSFile = fopen(MSFileName.c_str(),"w");
  if (!MSFile) {
    cerr << "Could not open the mesh size report file " << MSFileName << endl;
    throw;
  }

  fprintf(MSFile,"iter\t#nodes\t#edges\t#faces\t#regions\n");
  fprintf(MSFile,"init\t%d\t%d\t%d\t%d\n",
          M_numVertices(mesh),M_numEdges(mesh),M_numFaces(mesh),M_numRegions(mesh));
  fclose(MSFile);

  // ------------------------------------------------
  // setup debug tools
  // ------------------------------------------------

  int debugLevel = parameters.debugLevel;
  ma->setDebugLevel(debugLevel);

  // --- cross checks between two executions ---
  bool openJournal = parameters.openJournal;
  if ( openJournal ) {
    
    // Will hold a list of operations and checks in this execution
    ma->openJournal();

    // Will compare the list to another one (and abort if different)
    string refJournalName = parameters.referenceJournal;
    if ( !refJournalName.empty() ) {
      std::cout<<"Warning: comparing execution with the journal \'"
               << refJournalName << "\'\n";
      ma->setReferenceJournal(refJournalName);
    }
  }

  // --- slivers output ---
  bool sliverReports = parameters.sliverReports;
  if ( sliverReports ) {
    ma->enableSliverReports();
  }
  bool sliverOps = parameters.testSliverOperators;
  ma->testSliverOperators(sliverOps);

  displayMemoryUsage("Monitoring and debug tools set");

  // ------------------------------------------------
  // TEST MOBILE OBJECTS
  // ------------------------------------------------
  
  if ( !strcmp(parameters.task.c_str(),"MobileObject") ) {

    // initial mesh check
    if ( debugLevel >= 2 && !(ma->checkTheMesh()) )  ma->abort(__LINE__,__FILE__);

#ifndef PARALLEL
    ma->storeInitialCoordinates();

    displayMemoryUsage("Initial coordinates stored");

    // ---- build the mobile objects ----
    cout << "Building the mobile objects...\n";
    double cpu_mobj_0 = CPUTime();
    mobileObjectSet* objs = new mobileObjectSet();
    for (unsigned int iObj = 0; iObj < parameters.objects.size(); iObj++) {
      mobileObject* mob = new mobileObject(mesh);
      setMobileObject(mob,mesh,parameters.objects[iObj]);
      objs->insert(mob);
    }
    ma->registerObjects(objs);
    double cpu_mobj_tot = CPUTime() - cpu_mobj_0;
    cout << "Built the mobile objects in "<<cpu_mobj_tot<<" seconds\n";
    displayMemoryUsage("Mobile objects built");

// #warning "Size field printed to debug"
//     ma->writePos("sizeFieldInitMean.pos", OD_SIZEFIELD_MEAN);
//     ma->writePos("sizeFieldInitMin.pos",  OD_SIZEFIELD_MIN);
//     ma->writePos("sizeFieldInitMax.pos",  OD_SIZEFIELD_MAX);
//     ma->writePos("sizeFieldInit0.pos",OD_ANISO_SF_AXIS0);
//     ma->writePos("sizeFieldInit1.pos",OD_ANISO_SF_AXIS1);
//     ma->writePos("sizeFieldInit2.pos",OD_ANISO_SF_AXIS2);
// #warning "output curvature"
//     string nameCurv = "curvatureInit";
//     ma->writeVolumicCurvature(nameCurv);
// #warning "output distance"
//     string nameDist = "distanceInit";
//     ma->writeDistanceToWalls(nameDist);

    // ---- get node motion parameters ----
    string nodeMotionType = parameters.nodeMotionType;
    double elasticChi = parameters.elasticStiffnessAlteration;
    bool   elasticCavityMesh = parameters.elasticIsMeshTheCavity;
    int    elasticCavityThick = parameters.elasticCavityThickness;
#else
    string nodeMotionType = "None";
    double elasticChi;
    bool   elasticCavityMesh;
    int    elasticCavityThick;
#endif

    // ---- get time parameters ----
    double dtGlob = parameters.timeStep;
    double finalTime = parameters.finalTime;
    int maxTimeSteps = parameters.maxNbTimeSteps;



// #warning "hack propeller: to clean the CAD"
//     int cnt = 0;
//     EIter ei = M_edgeIter(mesh);
//     while ( pEdge pe = EIter_next(ei) ) {
//       if ( E_whatInType(pe) == 1 ) {
//         if ( GEN_tag(E_whatIn(pe)) == 29 || 
//              GEN_tag(E_whatIn(pe)) == 31 || 
//              GEN_tag(E_whatIn(pe)) ==  37 || 
//              GEN_tag(E_whatIn(pe)) ==  51 || 
//              GEN_tag(E_whatIn(pe)) ==  43 || 
//              GEN_tag(E_whatIn(pe)) ==  49 || 
//              GEN_tag(E_whatIn(pe)) ==  57 || 
//              GEN_tag(E_whatIn(pe)) ==  59    )
//           {
//             ma->collapseEdgeBrute(pe);
//             cnt++;
//           }
//       }
//     } 
//     EIter_delete(ei);
//     printf("Brute collapsed %d edges\n",cnt);







    // ---- the loop ----
    double cpu_t0 = CPUTime();
    double t = 0., dt_iter = 0.;
    int iter = 0;
    double cpu_move = 0., cpu_adapt = 0., cpu_output = 0.;
    while (iter <= maxTimeSteps && t < finalTime)
      {
        double cpu_t1 = CPUTime();
        if ( iter > 0 )  // a first step is performed with the adaptation before any motion occurs
          {
            // --- find time step for this iteration ---
            double t_left_tot = std::max(finalTime - t, 0.0);
            dt_iter = std::min(dtGlob, t_left_tot);
          
            // --- move nodes (dummy move and Laplace smoothing) ---
            if ( !strcmp(nodeMotionType.c_str(),"ForceBoundaries") ) {

              std::cout << "Warning: forcing boundaries is followed by Laplace smoothing\n";
              throw;

              double part;
              if (!ma->partlyMoveObjects(t, dt_iter, &part)) {
                cout <<"Could not perform a forced move\n";
                ma->abort(__LINE__,__FILE__);
              }
              dt_iter *= part;

              if ( debugLevel >= 1 && !(ma->checkTheMesh()) )  ma->abort(__LINE__,__FILE__);
              ma->LaplaceSmoothing();
              if ( debugLevel >= 1 && !(ma->checkTheMesh()) )  ma->abort(__LINE__,__FILE__);
            }

            t += dt_iter;  ma->setTime(t);
            printf("\nGlobal time step %d, t = %f, dt = %f\n\n",iter,t,dt_iter);

            // --- move nodes (elasticity analogy) ---
            if ( !strcmp(nodeMotionType.c_str(),"ElasticAnalogy") ) {
              ma->moveObjectsAndReposition(t, dt_iter, 
                                           parameters.elasticSubAdapt,
                                           parameters.elasticSubAdaptQualThr,
                                           elasticChi,
                                           elasticCavityMesh, 
                                           elasticCavityThick);
              displayMemoryUsage("Nodes repositionned");
            }
          }
        double cpu_t2 = CPUTime(); 
        cpu_move += ( cpu_t2 - cpu_t1 );

        // --- adapt the mesh ---
        if ( parameters.adapt || iter == 0 ) {
          if ( parameters.adaptAlways || iter == 0  ) {
            ma->run();
            displayMemoryUsage("After adaptation");
          }
          else {
            double worstQ, meanQ;
            ma->getStatistics(&meanQ,&worstQ);
            if ( worstQ < parameters.adaptQualityThr ) {
              ma->run();
              displayMemoryUsage("After adaptation");
            }
          }
        }
        double cpu_t3 = CPUTime();
        cpu_adapt += ( cpu_t3 - cpu_t2 );

        // --- outputs ---

        {  // full mesh
          if ( ( outputFrequency > 0 ) && 
               ( (iter%outputFrequency) == 0 ) ) {
            writeSolution(ma,iter,outputType);
          }
          displayMemoryUsage("After outputs");
        }

        {  // modifications monitor
          int nSplit, nColl, nSwap;
          double cpu_spl, cpu_coll, cpu_swap, cpu_sliv;
          ma->getModificationsInfo(&nSplit, &nColl, &nSwap,
                                   &cpu_spl, &cpu_coll, &cpu_swap, &cpu_sliv);
          FILE* sf = fopen(statFileName.c_str(),"a");
          fprintf(sf,"%d\t%d\t%d\t%d\t%e\t%e\t%e\t%e\n", iter, nSplit, nColl, nSwap,
                  cpu_spl, cpu_coll, cpu_swap, cpu_sliv);
          fclose(sf);
        }

        {  // quality monitor
          ma->getStatistics(&meanShape,&worstShape);
          FILE* qf = fopen(qualityFileName.c_str(),"a");
          fprintf(qf,"%d\t%e\t%e\n",iter,meanShape,worstShape);
          fclose(qf);
        }

        {  // cpu monitor
          double cpu_t4 = CPUTime();
          cpu_output += ( cpu_t4 - cpu_t3 );
      
          FILE* cpuF = fopen(cpuFileName.c_str(),"a");
          fprintf(cpuF,"%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                  iter, cpu_t2 - cpu_t1, cpu_t3 - cpu_t2, cpu_t4 - cpu_t3,
                  cpu_move, cpu_adapt, cpu_output, cpu_t4 - cpu_t0 );
          fclose(cpuF);
        }

        {  // mesh size monitor
          FILE* msf = fopen(MSFileName.c_str(),"a");
          fprintf(msf,"%d\t%d\t%d\t%d\t%d\n",
                  iter,M_numVertices(mesh),M_numEdges(mesh),M_numFaces(mesh),M_numRegions(mesh));
          fclose(msf);
        }

        { // slivers stats
          if ( (iter%10) == 0 ) {
            ofstream sliverOut( (outputPrefix + "slivers_tmp").c_str() );
            ma->printSliverRegionStatistics(sliverOut);
          }
        }

        { // operations stats
          if ( (iter%10) == 0 ) {
            ofstream operationsOut( (outputPrefix + "statistics_tmp").c_str() );
            ma->printStatistics(operationsOut);
          }
        }

        { // journal
          if ( (iter%10) == 0 ) {
            ofstream journalOut( (outputPrefix + "journal_tmp").c_str() );
            ma->flushJournal(journalOut);
          }
        }

        displayMemoryUsage("After monitors");
        iter++;
      }

#ifndef PARALLEL
    if (objs) {
      deleteObjects(objs);
      delete objs;
    }
    displayMemoryUsage("Objects deleted");

    ma->removeStoredCoordinates();
    displayMemoryUsage("Stored coordinates deleted");
#endif
  }

  // ------------------------------------------------
  // TEST REMOVE REGION
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"RemoveRegion") )
    {
      for (int i=0; i<1; i++) {
        int numOp = 0;
        int oriNumFaces = M_numFaces(mesh);
        int count = 0;
        FIter fi = M_faceIter(mesh);
        pFace pf;
        while ( ( pf = FIter_next(fi) ) && ( count < oriNumFaces ) ) {
          count++;
          if (F_numRegions(pf) != 1) continue;
          else {
            pRegion pr = F_region(pf,0);
            if (!pr) pr = F_region(pf,1);
            if ( ma->removeRegion(pr) ) numOp++;
          }
        }
        FIter_delete(fi);
        cout << "Num region removed: "<<numOp<<endl;
      }
    }

  // ------------------------------------------------
  // test size fields intersection or smoothing
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"SizeField") )
    {
      ma->writePos("testSF.pos",OD_SIZEFIELD_MEAN);
      ma->writePos("aniso0.pos",OD_ANISO_SF_AXIS0);
      ma->writePos("aniso1.pos",OD_ANISO_SF_AXIS1);
      ma->writePos("aniso2.pos",OD_ANISO_SF_AXIS2);

      PWLSField* testSF = new PWLSField(mesh);
      string h = "2.*x+0.1";
      AnalyticalSField* testASF = new AnalyticalSField(h);

      testSF->intersect(testASF);
//       testSF->smooth(0.5);
  
      MeshAdapter* testMA = new MeshAdapter(mesh,testSF);
  
      testMA->writePos("sizeField.pos",OD_SIZEFIELD_MEAN);
       //      testMA->writePos("sizeFieldSmoothed.pos",OD_SIZEFIELD_MEAN);

      exit(0);
    }

  // ------------------------------------------------
  // TEST REMOVE SLIVER FACES
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"RemoveSliverFaces") )
    {
      for (int i=0; i<10; i++) {
        ma->removeSlivers();
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__);
      }
    }

  // ------------------------------------------------
  // TEST LAPLACE SMOOTHING
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"LaplaceSmoothing") )
    {
      for (int i=0; i<1; i++) {
        ma->LaplaceSmoothing();
        cout << "Smoothing "<<i<<" operated\n";
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__);
      }
    }

  // ------------------------------------------------
  // TEST GEOMETRIC CONSTRAINTS
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"GeoConstraints") )
    {
      ma->setConstraint(2,5);
      ma->removeConstraint(2,5);
    }

  // ------------------------------------------------
  // TEST CONSTRAINTS
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"TopoConstraints") )
    {
      FIter fi = M_faceIter(mesh);
      while ( pFace pf = FIter_next(fi) ) {
        int gtag  = GEN_tag ( EN_whatIn ((pEntity)pf) );
        int gtype = GEN_type( EN_whatIn ((pEntity)pf) );
        if ( gtype == 2 && ( gtag == 5 || gtag == 22 || gtag == 14 || gtag == 27 ) ) ma->setConstraint((pEntity)pf);
      } 
      FIter_delete(fi);

      VIter vi = M_vertexIter(mesh);
      pVertex pv;
      while ( ( pv = VIter_next(vi) ) ) {
        // int gtag  = GEN_tag ( EN_whatIn ((pEntity)pv) );
        int gtype = GEN_type( EN_whatIn ((pEntity)pv) );
        if ( gtype == 0 ) ma->setConstraint((pEntity)pv);
      } 
      VIter_delete(vi);
    }

  // ------------------------------------------------
  // TEST FACE SWAP OPERATOR
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"FaceSwap") )
    {
      for (int i=0; i<5; i++) {
        int count=0, numOp=0;
        int oriNumFaces = M_numFaces(mesh);
        FIter fi = M_faceIter(mesh);
        pFace pf;
        while ( ( pf = FIter_next(fi) ) && ( count < oriNumFaces ) ) {
          int res = ma->swapFace(pf);
          if ( res ) numOp++;
          count++;
        }
        FIter_delete(fi);
        cout << "Num face swaps: "<<numOp<<endl;
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__);

        int nESwaps = ma->optimiseElementShape();
        cout<< "Edge swaps applied: "<<nESwaps<<endl;
      }
    }

  // ------------------------------------------------
  // MAKE A BAD MESH
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"UglyMesh") )
    {
      ma->uglyTheMesh(0.2,10);
      ma->writeMsh("uglyMesh.msh");
    }

  // ------------------------------------------------
  // TEST REMOVE SLIVER REGIONS
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"RemoveSliverRegions") )
    {
      ma->setSliverQuality(0.1);
      for (int i=0; i<2; i++) {
        ma->removeSlivers();
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__);
      }
    }
  
  // ------------------------------------------------
  // TEST DOUBLE-SPLIT-COLLAPSE OPERATOR
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"DESC") )
    {
      for (int i=0; i<20; i++) {
        int numOp = 0;
        int oriNumRgn = M_numRegions(mesh);
        int count = 0;
        RIter ri = M_regionIter(mesh);
        pRegion pr;
        while ( ( pr = RIter_next(ri) ) && ( count < oriNumRgn ) ) {
          pEdge pe1 = R_edge(pr,0);
          pEdge pe2 = R_gtOppEdg(pr, pe1);

          if ( ma->DSplitCollapseEdge(pr,pe1,pe2) )
            numOp++;
          count++;
        }
        RIter_delete(ri);
        cout << "Num double-split-collapses: "<<numOp<<endl;

        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__);
        //     ma->printStatistics(std::cout);

        //     int* nbBefore = new int(0);
        //     int* nbAfter = new int(0);
        //     ma->removeSliverRegions(nbBefore,nbAfter);
        //     cout <<"Nb slivers before->after realignement:\t"<<*nbBefore<<" -> "<<*nbAfter<<endl;
      }
    }

  // ------------------------------------------------
  // TEST EDGE SWAP
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"EdgeSwap") )
    {
      if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__); 
      for (int i=0; i<30; i++) {
        int count = 0;
        int numOp = 0;
        int numE = M_numEdges(mesh);
        EIter ei = M_edgeIter(mesh);
        while ( pEdge pe = EIter_next(ei) ) {
          count++;
          if ( ma->swapEdge(pe) ) numOp++;
          if ( count > numE ) break;
        } 
        EIter_delete(ei);
        cout << "Num edge swaps ("<<i<<"): " << numOp << endl;
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__); 
        writeSolution(ma,i,outputType);
      }
    }

  // ------------------------------------------------
  // TEST GEOMATCHER
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"GeoMatcher") )
    {
      for (int i=0; i<8; i++) {
        int count = 0;
        int numOp = 0;
        int numE = M_numEdges(mesh);
        EIter ei = M_edgeIter(mesh);
        pEdge pe;
        while ( ( pe = EIter_next(ei) ) ) {
          count++;
          if ( ma->splitEdge(pe) ) {
            ma->snapVertices();
            numOp++;
          }
          if ( count > numE ) break;
        } 
        EIter_delete(ei);
        cout << "Num edge splits ("<<i<<"): " << numOp << endl;
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__); 
      }
    }

  // ------------------------------------------------
  // TEST EDGE SPLIT
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"EdgeSplit") )
    {
      if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__); 
      for (int i=0; i<5; i++) {

        ma->splitEveryEdgeOnce();

//         int count = 0;
//         int numOp = 0;
//         int numE = M_numEdges(mesh);
//         EIter ei = M_edgeIter(mesh);
//         pEdge pe;
//         while ( ( pe = EIter_next(ei) ) ) {
//           count++;
//           if ( E_whatInType(pe) <= 2 && ma->splitEdge(pe) ) { numOp++; }
//           if ( count > numE ) break;
//         } 
//         EIter_delete(ei);
//         cout << "Num edge splits ("<<i<<"): " << numOp << endl;
//         if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__); 

//         // --- Geometry tracking ---
//         cout << "Snapping boundary nodes" << endl;
//         if (numOp) ma->snapVertices();

        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__);
        writeSolution(ma,i,outputType);
      }
    }

  // ------------------------------------------------
  // TEST EDGE COLLAPSE
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"EdgeCollapse") )
    {
      //      GM_writeGEO(model,"testgeo.geo");
      //      ma->writeMsh("MeshInit.msh");
      for (int i=0; i<10; i++) {
        double t0 = CPUTime();
        int count = 0;
        int nbEdgesInit = M_numEdges(mesh);
        int numOp = 0;
        EIter ei = M_edgeIter(mesh);
        pEdge pe;
        while ( ( pe = EIter_next(ei) ) ) {
          count++;
          if ( ma->collapseEdge(pe) ) { numOp++; }
          if ( count > nbEdgesInit ) break;
        } 
        EIter_delete(ei);
        double dt = CPUTime() - t0;
        cout << "Num edge collapses ("<<i<<"): " << numOp << " in "<<dt<<" seconds\n";
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__); 
//         ma->printStatistics(std::cout);
        writeSolution(ma,i,outputType);
      }
    }

  // ------------------------------------------------
  // TEST EDGE SPLIT AND COLLAPSE
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"EdgeSplitCollapse") )
    {
      for (int i=0; i<20; i++) {
        
        // --- collapse as many edges as possible ---
        int totColl = 0;
        while (1) {
          int numClp = 0;
//           while (1)
//             {
//               pEdge pe;
//               EIter ei = M_edgeIter(mesh);
//               while ( ( pe = EIter_next(ei) ) ) {
//                 if ( ma->collapseEdge(pe) ) { numClp++; break; }
//               }
//               EIter_delete(ei);
//               if ( !pe ) break;
//               if ( !checkMesh(mesh,CHECK_GEOM_COMPATIBILITY,1,std::cout) ) ma->abort(__LINE__,__FILE__); 
//             }
          EIter ei = M_edgeIter(mesh);
          pEdge pe;
          while ( ( pe = EIter_next(ei) ) ) {
            if ( ma->collapseEdge(pe) ) numClp++;
          }
          EIter_delete(ei);
          totColl += numClp;
          cout << "Num edge collapses: "<< numClp << endl;
          if ( ! ma->checkTheMesh() ) ma->abort(__LINE__,__FILE__); 
          if ( numClp == 0 ) break;
        }

        // --- perform many edge splits ---
        for (int iS=0; iS<2; iS++) {
          int numSpl = ma->splitEveryEdgeOnce();
          cout << "Num edge splits: " << numSpl << endl;
          if ( ! ma->checkTheMesh() ) ma->abort(__LINE__,__FILE__); 
        }

        // --- Geometry tracking ---
//         if (numSpl) ma->snapVertices();
//         if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__); 

        writeSolution(ma,i,outputType);
      }
    }

  // ------------------------------------------------
  // TEST EDGE SPLIT AND COLLAPSE 2
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"EdgeSplitCollapse2") )
    {
      for (int i=0; i<10; i++) {
        int count = 0;
        int oriNumEdges = M_numEdges(mesh);
        int numClp = 0, numSpl = 0;
        EIter ei = M_edgeIter(mesh);
        pEdge pe;
        while ( ( pe = EIter_next(ei) ) ) {
          count++;
          if ( E_whatInType(pe) <= 2 ) 
            {
              if ( (count%3) == 0 ) {
                if ( ma->splitEdge(pe) ) numSpl++;
              }
              else {
                if ( ma->collapseEdge(pe) ) numClp++;
              }
            }
          if ( count > oriNumEdges ) break;
        } 
        EIter_delete(ei);
        cout << "Num edge splits ("<<i<<"): " << numSpl << endl;
        cout << "Num edge collapses ("<<i<<"): " << numClp << endl;
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__); 

        // --- Geometry tracking ---
        if (numSpl) ma->snapVertices();
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__); 

        writeSolution(ma,i,outputType);
//         ma->printStatistics(std::cout);
      }
      
    }

  // ------------------------------------------------
  // TEST FACE COLLAPSE OPERATOR
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"FaceCollapse") )
    {
      for (int i=0; i<10; i++) {
        int numOp = 0;
        int oriNumFaces = M_numFaces(mesh);
        int count = 0;
        FIter fi = M_faceIter(mesh);
        pFace pf;
        while ( ( pf = FIter_next(fi) ) && ( count < oriNumFaces ) ) {
          for (int j=0; j<3; j++) {
            pEdge pe = F_edge(pf,j);
            if ( ma->collapseFace(pf,pe) ) {
              numOp++;
              break;
            }
          }
          count++;
        }
        FIter_delete(fi);
        cout << "Num face collapse: "<<numOp<<endl;
        if ( !checkMesh(mesh,CHECK_ALL,1) )  ma->abort(__LINE__,__FILE__);
//         ma->printStatistics(std::cout);
      }
    }

  // ------------------------------------------------
  // TEST OPTIMISE EDGE LENGTH
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"OptimizeLength") )
    {
      cout<<"Split collapses applied: "<<ma->optimiseEdgeLength()<<endl;
    }

  // ------------------------------------------------
  // TEST OPTIMISE SHAPE
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"OptimizeShape") )
    {
      int nSwaps = ma->optimiseElementShape();
      cout<< "Swaps applied: "<<nSwaps<<endl;
    }

  // ------------------------------------------------
  // TEST METRIC
  // ------------------------------------------------

//   else if ( !strcmp(parameters.task.c_str(),"Metric") )
//     {
//       double h[3] = { 2., 3., 10. };
//       double e[3][3] = {
//         {sqrt(2.)/2., sqrt(2.)/2., 0.},
//         {0., 0., 1.},
//         {sqrt(2.)/2., -sqrt(2.)/2., 0.}
//       };

//       MAdMetric M = MAdMetric(h[0],h[1],h[2],e[0],e[1],e[2]);

//       doubleMatrix V = doubleMatrix(3,3);
//       doubleVector S = doubleVector(3);
//       M.eig(V,S,true);

//       printVec(h,"h");
//       printMat(e,"e");
//       V.print("V");
//       S.print("S");
//     }

  // ------------------------------------------------
  // TEST CURVATURES
  // ------------------------------------------------

  else if ( !strcmp(parameters.task.c_str(),"Curvatures") )
    {
      ma->writePos("CurvDiv.pos",OD_CURVATURE_DIV);
      ma->writePos("CurvMax.pos",OD_CURVATURE_MAX);
      ma->writePos("CurvMin.pos",OD_CURVATURE_MIN);
      ma->writePos("CurvMaxVec.pos",OD_CURVATURE_MAX_VEC);
      ma->writePos("CurvMinVec.pos",OD_CURVATURE_MIN_VEC);
      
//       FIter fi = M_faceIter(mesh);
//       pFace pf;
//       while ( ( pf = FIter_next(fi) ) )
//         {
//           if ( F_whatInType(pf) != 2 ) continue;
          
//           pGEntity pge = F_whatIn(pf);
//           {
//             double u[2] = { 0., 0. };
//             double xyz[3];
//             GF_xyz( (pGFace) pge, u[0], u[1], xyz);
//             printf("u: %f, v: %f\n",u[0],u[1]);
//             printVec(xyz,"xyz");
//           }
          
//           double u[2] = { 0.5, 0.5 };
//           double xyz[3];
//           GF_xyz( (pGFace) pge, u[0], u[1], xyz);
//           printf("u: %f, v: %f\n",u[0],u[1]);
//           printVec(xyz,"xyz");
          
//           double divCurv = GF_curvature( (pGFace)pge, u);
//           double dirMax[3], dirMin[3], maxCurv, minCurv;
//           double mCurv = GF_curvatures( (pGFace)pge, u,
//                                         dirMax, dirMin, &maxCurv, &minCurv);
          
//           printVec(dirMax,"dirMax");
//           printVec(dirMin,"dirMin");
//           printf("Max curv: %f, min curv: %f, div curv: %f\n",maxCurv,minCurv,divCurv);
          
//           double tmp[3];
//           tmp[0] = dotProd(xyz,dirMax);
//           tmp[1] = dotProd(xyz,dirMin);
//           tmp[2] = dotProd(dirMin,dirMax);
//           printVec(tmp,"sould be 0 vector");
          
//           break;
//         }
//       FIter_delete(fi);
    }

  // ------------------------------------------------
  // MISCELLANEOUS
  // ------------------------------------------------
  
  else if ( !strcmp(parameters.task.c_str(),"Miscellaneous") )
    {
      int count = 0;
      EIter ei = M_edgeIter(mesh);
      pEdge pe;
      while ( ( pe = EIter_next(ei) ) ) {
        count++;
        if ( count == 69 ) //65
          {
            if ( !(ma->collapseEdge(pe)) ) {
              printf("not collapsed\n");
              throw;
            }
            break;
          }
      } 
      EIter_delete(ei);
    }

  // ------------------------------------------------
  else if ( !strcmp(parameters.task.c_str(),"None") )
    {
    }

  else {
    std::cerr << "Error: unknown task " << parameters.task.c_str() << "\n";
  }

  // ------------------------------------------------
  // Write final result
  // ------------------------------------------------

  ma->writeMsh("result.msh");
//   ma->writePos("result.pos",OD_MEANRATIO);
  
// #warning "Size field printed to debug"
//     ma->writePos("sizeFieldFinalMean.pos", OD_SIZEFIELD_MEAN);
//     ma->writePos("sizeFieldFinalMin.pos",  OD_SIZEFIELD_MIN);
//     ma->writePos("sizeFieldFinalMax.pos",  OD_SIZEFIELD_MAX);
//     ma->writePos("sizeFieldFinal0.pos",OD_ANISO_SF_AXIS0);
//     ma->writePos("sizeFieldFinal1.pos",OD_ANISO_SF_AXIS1);
//     ma->writePos("sizeFieldFinal2.pos",OD_ANISO_SF_AXIS2);
// #warning "output curvature"
//     ma->writeVolumicCurvature("curvatureFinal");
// #warning "output distance"
//     ma->writeDistanceToWalls("distanceFinal");

  ma->printStatistics(cout);
  ofstream statOut((outputPrefix + "statistics").c_str());
  ma->printStatistics(statOut);

  ofstream sliverOut((outputPrefix + "slivers").c_str());
  ma->printSliverRegionStatistics(sliverOut);

  ofstream journalOut((outputPrefix + "journal").c_str());
  ma->flushJournal(journalOut);

  // ------------------------------------------------
  // Clean up
  // ------------------------------------------------

  if (ma) delete ma;
  set<pSField>::iterator sfIter = sizeFields.begin();
  set<pSField>::iterator sfLast = sizeFields.end();
  for (; sfIter != sfLast; sfIter++) {
    delete (*sfIter);
  }
  displayMemoryUsage("Adapter and size fields deleted");
  
  if (mesh)  M_delete(mesh);
  if (model) GM_delete(model);
  displayMemoryUsage("Mesh and model deleted");

  MAdLibFinalize();

  printf("Total execution time: %f seconds\n", CPUTime()-main_t0);
}

// ----------------------------------------------------------------------
