/* a semple converter between gmsh and excudus file format */

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include "Gmsh.h"
#include "GModel.h"
#include "GVertex.h"
#include "MElement.h"
#include "exodusII.h"

typedef std::set<int> tri;
typedef std::vector<tri> triVec;

void errExit(int errCode, std::string msg)
{
  if (errCode) {
    std::cout<< "Error!" << msg 
	     << " with error code "
	     << errCode << std::endl;
    exit(-1);
  }
  return;
}

bool isIn(int qry, std::vector<int>& vec)
{
  bool found = false;
  for (auto it=vec.begin(); it!=vec.end(); it++)
    if (*it == qry)
      return(true);
  return(false);
}

bool isIn(tri& qry, triVec& vec)
{
  for (auto it=vec.begin(); it!=vec.end(); it++)
    if (*it == qry)
      return(true);
  return(false);
}

int main(int argc, char* argv[])
{
  std::cout << "WARNING! This converter is deprecated and will be removed in future releases.\n";
  std::cout << "WaRNING! We recommend using nemosysRun's conversion service.\n";
  std::cout << "WaRNING! This converter only works for tetrahedral elements.\n";
  if (argc < 5){
    std::cout << "Gmsh to Exodus II converter v1.0\n";
    std::cout << "Converter between gmsh and exodusII v1.0\n";
    std::cout << "Usage: \n";
    std::cout << argv[0] << " gmshProjectileFile outputExoFile nTargetFile  gmshTargetFile1 ...\n\n";
    return(0);
  }
  int exErr;
  int exo_fid;
  int nTrg;

  std::string gmshPrjFName, exoFName;
  gmshPrjFName = argv[1];
  exoFName = argv[2];
  nTrg = std::stoi(argv[3]);
  std::vector<std::string> gmshTrgFName;
  for (int iTrg=0; iTrg<nTrg; iTrg++)
    gmshTrgFName.push_back(std::string(argv[4+iTrg]));
 
  // read gmsh file
  std::cout << "Reading projectile and target mesh files.\n";
  GModel *prjGM, *trgGM;
  prjGM= new GModel("projectile");
  prjGM->readMSH(gmshPrjFName);
  std::cout << "Number of nodes\n";
  std::cout << "  projectile = " << prjGM->getMaxVertexNumber() << std::endl;

  ///////////////////////////////////////////////////////////////////////////
  // find nodal coordinate
  std::cout << "Reading nodal coordinates.\n";
  std::vector<double> x,y,z;
  for (int iVrt=1; iVrt <= prjGM->getMaxVertexNumber(); iVrt++)
  {
    MVertex* vrt = prjGM->getMeshVertexByTag(iVrt);
    x.push_back(vrt->x());
    y.push_back(vrt->y());
    z.push_back(vrt->z());
  }

  ///////////////////////////////////////////////////////////////////////////
  // getting connectivities 
  // we assume just tet elements
  std::cout << "Reading element connectivity.\n";
  std::vector<int> elmConnPrjGMSH;
  std::vector<int> elmConnPrjTriGMSH;
  triVec triConnPrj;
  int nTetPrj = 0;
  int nTriPrj = 0;
  for (int iElm=1; iElm<=prjGM->getNumMeshElements(); iElm++)
  {
    std::vector<MVertex*> verts;
    MElement* elm = prjGM->getMeshElementByTag(iElm);
    if (elm->getType() == 5)
    {
      nTetPrj++;
      elm->getVertices(verts);
      for (auto iv=verts.begin(); iv!=verts.end(); iv++)
      {
	elmConnPrjGMSH.push_back((*iv)->getNum());
	//std::cout << (*iv)->getNum() << std::endl;
      }
    } else if (elm->getType() == 3) {
      nTriPrj++;
      elm->getVertices(verts);
      for (auto iv=verts.begin(); iv!=verts.end(); iv++)
      {
	elmConnPrjTriGMSH.push_back((*iv)->getNum());
	//std::cout << (*iv)->getNum() << std::endl;
      }
      tri tmp = {verts[0]->getNum(), verts[1]->getNum(), verts[2]->getNum()};
      triConnPrj.push_back(tmp);
    }
  }
  
  // working on targets
  int totNdeTrg = 0;
  std::vector<int> nNdeTrg;
  nNdeTrg.resize(nTrg,0);
  int totTetTrg = 0;
  std::vector<int> nTetTrg;
  nTetTrg.resize(nTrg,0);
  int nTriTrg = 0;
  std::vector<std::vector<int> > elmConnTrgGMSH;
  triVec triConnTrg;
  int totNNdeSoFar = prjGM->getMaxVertexNumber();
  for (int iTrg=0; iTrg<nTrg; iTrg++)
  {
    trgGM= new GModel(gmshTrgFName[iTrg].c_str());
    trgGM->readMSH(gmshTrgFName[iTrg].c_str());
    std::cout << "Number of nodes " << gmshTrgFName[iTrg];
    std::cout <<  " " << trgGM->getMaxVertexNumber() << std::endl;
    nNdeTrg[iTrg] = trgGM->getMaxVertexNumber();
 
    /* 
    ///////////////////////////////////////////////////////////////////////////
    // get physical group dim=2, num=1 and find all vertices residing in it
    std::vector<MVertex*> prjBndryMVrtxObj;
    prjGM->getMeshVerticesForPhysicalGroup(2, 1, prjBndryMVrtxObj);
    std::vector<int> prjBndryVrtxId;
    for (auto it=prjBndryMVrtxObj.begin(); it!=prjBndryMVrtxObj.end(); it++)
    {
      prjBndryVrtxId.push_back((*it)->getNum());
    }
    */

    ///////////////////////////////////////////////////////////////////////////
    // attach nodal coordinate
    std::cout << "Reading nodal coordinates.\n";
    for (int iVrt=1; iVrt <= trgGM->getMaxVertexNumber(); iVrt++)
    {
      totNdeTrg++;
      MVertex* vrt = trgGM->getMeshVertexByTag(iVrt);
      x.push_back(vrt->x());
      y.push_back(vrt->y());
      z.push_back(vrt->z());
    }

    ///////////////////////////////////////////////////////////////////////////
    // getting connectivities 
    // we assume just tet elements
    std::vector<int> tmpElmConn;
    std::cout << "Reading element connectivity.\n";
    for (int iElm=1; iElm<=trgGM->getNumMeshElements(); iElm++)
    {
      std::vector<MVertex*> verts;
      MElement* elm = trgGM->getMeshElementByTag(iElm);
      if (elm->getType() == 5)
      {
	nTetTrg[iTrg]++;
        totTetTrg++;
	elm->getVertices(verts);
	for (auto iv=verts.begin(); iv!=verts.end(); iv++)
	{
	  tmpElmConn.push_back((*iv)->getNum() + totNNdeSoFar);
	  //std::cout << (*iv)->getNum() << std::endl;
	}
      } else if (elm->getType() == 3) {
	nTriTrg++;
	elm->getVertices(verts);
	tri tmp = {verts[0]->getNum(), verts[1]->getNum(), verts[2]->getNum()};
	triConnTrg.push_back(tmp);
      }
    }
    std::cout << "Number of tetra elements ";
    std::cout << nTetTrg[iTrg] << std::endl;
    elmConnTrgGMSH.push_back(tmpElmConn);
    totNNdeSoFar += trgGM->getMaxVertexNumber();
  
    /*
    ///////////////////////////////////////////////////////////////////////////
    // getting boundary face id element id pairs
    // we assume just tet elements
    int fG[4][3] = {{0,1,3}, {1,2,3}, {0,2,3}, {0,1,2}};
    std::vector<int> bndryElmIdxPrj;
    std::vector<int> bndryFaceIdxPrj;
    int cntTet = 0;
    int nSideTet = 0;
    for (int iElm=1; iElm<=prjGM->getNumMeshElements(); iElm++)
    {
      MElement* elm = prjGM->getMeshElementByTag(iElm);
      if (elm->getType() == 5)
      {
	cntTet++;
	std::vector<MVertex*> vs;
	elm->getVertices(vs);
	for (int iF=0; iF<4; iF++)
	{
	  tri tmp = {vs[fG[iF][0]]->getNum(), vs[fG[iF][1]]->getNum(), vs[fG[iF][2]]->getNum()};
	  if (isIn(tmp, triConnPrj))
	  {
	    nSideTet++;
	    bndryElmIdxPrj.push_back(cntTet);
	    bndryFaceIdxPrj.push_back(iF+1);
	  }
		      
	}
      }
    }
    std::cout << "Number of projectile side tets = " << nSideTet << std::endl;
    // Target
    std::vector<int> bndryElmIdxTrg;
    std::vector<int> bndryFaceIdxTrg;
    cntTet = 0;
    nSideTet = 0;
    for (int iElm=1; iElm<=trgGM->getNumMeshElements(); iElm++)
    {
      MElement* elm = trgGM->getMeshElementByTag(iElm);
      if (elm->getType() == 5)
      {
	cntTet++;
	std::vector<MVertex*> vs;
	elm->getVertices(vs);
	for (int iF=0; iF<4; iF++)
	{
	  tri tmp = {vs[fG[iF][0]]->getNum(), vs[fG[iF][1]]->getNum(), vs[fG[iF][2]]->getNum()};
	  if (isIn(tmp, triConnTrg))
	  {
	    nSideTet++;
	    bndryElmIdxTrg.push_back(cntTet);
	    bndryFaceIdxTrg.push_back(iF+1);
	  }
		      
	}
      }
    }
    std::cout << "Number of target side tets = " << nSideTet << std::endl;
    */
  }

  ///////////////////////////////////////////////////////////////////////////
  std::cout << "Creating exodus II output file.\n";
  int comp_ws = 8;
  int io_ws = 8;
  // create exo file
  exo_fid = ex_create(exoFName.c_str(), EX_CLOBBER, &comp_ws, &io_ws);
  // initialize the file
  exErr = ex_put_init(exo_fid, "gmsh", 3, 
              x.size(), 
              prjGM->getNumMeshElements() + totTetTrg, 
              1+nTrg, 2, 0);
  errExit(exErr, "Error in initialization"); 

  ///////////////////////////////////////////////////////////////////////////
  // write nodal coordinates
  std::cout << "Writing nodal coordinates.\n";
  exErr = ex_put_coord(exo_fid, &x[0], &y[0], &z[0]);
  errExit(exErr, "Error in writing node coordinate"); 

  ///////////////////////////////////////////////////////////////////////////
  // write node set for the nodes of projectile
  std::cout << "Writing node set.\n";
  exErr = ex_put_node_set_param(exo_fid, 1, prjGM->getMaxVertexNumber(), 0);
  errExit(exErr, "Error in defining node set"); 
  std::vector<int> nodeSet;
  for (int iNde=1; iNde<=prjGM->getMaxVertexNumber(); iNde++)
    nodeSet.push_back(iNde);
  exErr = ex_put_node_set(exo_fid, 1, &nodeSet[0]);
  
  // write node set for the nodes of Target
  exErr = ex_put_node_set_param(exo_fid, 2, totNdeTrg, 0);
  errExit(exErr, "Error in defining node set"); 
  nodeSet.clear();
  for (int iNde=1; iNde<=totNdeTrg; iNde++)
    nodeSet.push_back(iNde + prjGM->getMaxVertexNumber());
  exErr = ex_put_node_set(exo_fid, 2, &nodeSet[0]);
  

  ///////////////////////////////////////////////////////////////////////////
  // write element block for projectile
  std::cout << "Writing element block for projectile.\n";
  exErr = ex_put_elem_block(exo_fid, 1, "TETRA", nTetPrj, 4, 1);
  errExit(exErr, "Error in writing element block");
  // write connectivity (node map)
  
  ///////////////////////////////////////////////////////////////////////////
  // performing additional tests on the projectile mesh 
  // to ensure compatability with EPIC requirements
  
  // what is the lowest and highest node index (should be 1 and prjMesh->getMaxVertexNumber() )
  std::cout << "Min nde indx prj = " << *min_element(elmConnPrjGMSH.begin(), elmConnPrjGMSH.end()) << "\n"
            << "Max nde indx prj = " << *max_element(elmConnPrjGMSH.begin(), elmConnPrjGMSH.end()) << "\n";
  // creating unique element id list and ensuring its size is equal to the number of projectile 
  // nodes, in other words there should not be any additional nodes in the projectile mesh that
  // are not referenced in the connectivity table of the tetrahedral elements
  std::set<int> unqPrjConn;
  for (auto it = elmConnPrjGMSH.begin(); it != elmConnPrjGMSH.end(); it++)
      unqPrjConn.insert(*it);
  if ( unqPrjConn.size() != prjGM->getMaxVertexNumber() )
  {
    std::cout << "Size uniqe node id table for the projectile = " << unqPrjConn.size() << std::endl;
    std::cerr << "The projectile mesh contains nodes that are not referenced by any tetrahedral element!\n";
    std::cerr << "Make sure projectile mesh only contains tetrahedral elements.\n";
    exit(-1);
  }

  exErr = ex_put_elem_conn(exo_fid, 1, &elmConnPrjGMSH[0]);
  errExit(exErr, "Error in writing element connectivities");
  // write element attribute
  std::vector<double> elmAttrib;
  elmAttrib.resize(prjGM->getNumMeshElements(), 1);
  exErr = ex_put_elem_attr(exo_fid, 1, &elmAttrib[0]);
  errExit(exErr, "Error in writing element attribute.");

  ///////////////////////////////////////////////////////////////////////////
  // write element block for target
  std::cout << "Writing element block for target(s).\n";
  for (int iTrg=0; iTrg<nTrg; iTrg++)
  {
    exErr = ex_put_elem_block(exo_fid, 2+iTrg, "TETRA", nTetTrg[iTrg], 4, 1);
    errExit(exErr, "Error in writing element block");
    // write connectivity (node map)
    exErr = ex_put_elem_conn(exo_fid, 2+iTrg, &elmConnTrgGMSH[iTrg][0]);
    errExit(exErr, "Error in writing element connectivities");
    // write element attribute
    elmAttrib.resize(nTetTrg[iTrg], 2);
    exErr = ex_put_elem_attr(exo_fid, 2+iTrg, &elmAttrib[0]);
    errExit(exErr, "Error in writing element attribute.");
  }

  /*
  ///////////////////////////////////////////////////////////////////////////
  // write TRI element block for projectile
  std::cout << "Writing TRI element block for projectile.\n";
  exErr = ex_put_elem_block(exo_fid, 3, "TRIANGLE", nTriPrj, 3, 1);
  errExit(exErr, "Error in writing element block");
  // write connectivity (node map)
  exErr = ex_put_elem_conn(exo_fid, 3, &elmConnPrjTriGMSH[0]);
  errExit(exErr, "Error in writing element connectivities");
  // write element attribute
  elmAttrib.resize(nTriPrj, 1);
  exErr = ex_put_elem_attr(exo_fid, 3, &elmAttrib[0]);
  errExit(exErr, "Error in writing element attribute.");

  ///////////////////////////////////////////////////////////////////////////
  // write side set for projectile
  exErr = ex_put_side_set_param(exo_fid, 1, bndryFaceIdxPrj.size(), 0);
  errExit(exErr, "Error in defining projectile side set.");
  exErr = ex_put_side_set(exo_fid, 1, &bndryElmIdxPrj[0], &bndryFaceIdxPrj[0]);
  errExit(exErr, "Error in writing projectile side set.");
  // write side set for target
  exErr = ex_put_side_set_param(exo_fid, 2, bndryFaceIdxTrg.size(), 0);
  errExit(exErr, "Error in defining target side set.");
  exErr = ex_put_side_set(exo_fid, 2, &bndryElmIdxTrg[0], &bndryFaceIdxTrg[0]);
  errExit(exErr, "Error in writing projectile side set.");
  */
  
  // ending
  ex_close(exo_fid);
  std::cout << "Application ended successfully\n";
  return(0);
}
