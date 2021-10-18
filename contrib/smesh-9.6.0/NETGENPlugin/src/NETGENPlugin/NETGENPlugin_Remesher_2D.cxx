// Copyright (C) 2007-2021  CEA/DEN, EDF R&D, OPEN CASCADE
//
// Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// File      : NETGENPlugin_Remesher_2D.cxx
// Created   : Thu Sep 21 16:48:46 2017
// Author    : Edward AGAPOV (eap)
//

#include "NETGENPlugin_Remesher_2D.hxx"

#include "NETGENPlugin_Mesher.hxx"
#include "NETGENPlugin_Hypothesis_2D.hxx"

#include <SMDS_SetIterator.hxx>
#include <SMESHDS_Group.hxx>
#include <SMESHDS_Mesh.hxx>
#include <SMESH_ControlsDef.hxx>
#include <SMESH_Gen.hxx>
#include <SMESH_MeshAlgos.hxx>
#include <SMESH_MesherHelper.hxx>
#include <SMESH_Group.hxx>
#include <SMESH_MeshEditor.hxx>
#include <SMESH_subMesh.hxx>

#include <Bnd_B3d.hxx>
#include <Precision.hxx>

#include <occgeom.hpp>
#include <meshing.hpp>
#include <stlgeom.hpp>
//#include <stltool.hpp>

#include <boost/container/flat_set.hpp>

using namespace nglib;

namespace netgen {

  NETGENPLUGIN_DLL_HEADER
  extern MeshingParameters mparam;

  NETGENPLUGIN_DLL_HEADER
  extern STLParameters stlparam;

  NETGENPLUGIN_DLL_HEADER
  extern STLDoctorParams stldoctor;
}
namespace nglib
{
  NETGENPLUGIN_DLL_HEADER
  extern netgen::Array<netgen::Point<3> > readedges;
}

namespace
{
  //=============================================================================
  /*!
   * \brief Fill holes in the mesh, since netgen can remesh only a closed shell mesh.
   *        At destruction, remove triangles filling the holes
   */
  class HoleFiller
  {
  public:
    HoleFiller( SMESH_Mesh& meshDS );
    ~HoleFiller();
    void AddHoleBordersAndEdges( Ng_STL_Geometry * ngStlGeo, bool toAddEdges );
    void KeepHole() { myHole.clear(); myCapElems.clear(); }
    void ClearCapElements() { myCapElems.clear(); }

  private:
    SMESHDS_Mesh*                          myMeshDS;
    std::vector< std::vector< gp_XYZ > >   myHole;      // initial border nodes
    std::vector< gp_XYZ >                  myInHolePos; // position inside each hole
    std::vector< const SMDS_MeshElement* > myCapElems;  // elements closing holes
  };

  //================================================================================
  /*!
   * \brief Fill holes in the mesh
   */
  //================================================================================

  HoleFiller::HoleFiller( SMESH_Mesh& theMesh ):
    myMeshDS( theMesh.GetMeshDS() )
  {
    SMESH_MeshEditor editor( &theMesh );
    {
      // // merge nodes
      // const double tol = Max( 0.1 * netgen::mparam.minh, Precision::Confusion() );
      // TIDSortedNodeSet allNodes;
      // SMESH_MeshEditor::TListOfListOfNodes equalNodes;
      // editor.FindCoincidentNodes( allNodes, tol, equalNodes, true );
      // editor.MergeNodes( equalNodes, /*noHoles=*/false );
    }

    // find holes
    SMESH_MeshAlgos::TFreeBorderVec holes;
    bool isManifold = true, isGoodOri = true;
    SMESH_MeshAlgos::FindFreeBorders( *myMeshDS, holes, /*closedOnly=*/true,
                                      &isManifold, &isGoodOri );

    if ( !isManifold )
    {
      // set bad faces into a compute error
      const char* text = "Non-manifold mesh. Only manifold mesh can be re-meshed";
      SMESH_BadInputElements* error =
        new SMESH_BadInputElements( myMeshDS, COMPERR_BAD_INPUT_MESH, text );
      SMESH::Controls::MultiConnection2D fun;
      fun.SetMesh( myMeshDS );
      SMDS_ElemIteratorPtr fIt = myMeshDS->elementsIterator( SMDSAbs_Face );
      while ( fIt->more() )
      {
        const SMDS_MeshElement* f = fIt->next();
        if ( fun.GetValue( f->GetID() ) > 2 )
          error->myBadElements.push_back( f );
      }
      theMesh.GetSubMesh( theMesh.GetShapeToMesh() )->GetComputeError().reset( error );
      throw SALOME_Exception( text );
    }

    // don't want to sew coincident borders
    if ( !holes.empty() )
    {
      // define tolerance
      double tol, len, sumLen = 0, minLen = 1e100;
      int     nbSeg = 0;
      for ( size_t i = 0; i < holes.size(); ++i )
      {
        nbSeg += holes[i].size();
        SMESH_NodeXYZ p1 = holes[i][0];
        for ( size_t iP = 1; iP < holes[i].size(); ++iP )
        {
          SMESH_NodeXYZ p2 = holes[i][iP];
          len = ( p1 - p2 ).Modulus();
          sumLen += len;
          minLen = Min( minLen, len );
          p1 = p2;
        }
      }
      double avgLen = sumLen / nbSeg;
      if ( minLen > 1e-5 * avgLen )
        tol = 0.1 * minLen; // minLen is not degenerate
      else
        tol = 0.1 * avgLen;

      SMESH_MeshAlgos::CoincidentFreeBorders freeBords;
      SMESH_MeshAlgos::FindCoincidentFreeBorders( *myMeshDS, tol, freeBords );
      if ( !freeBords._coincidentGroups.empty() )
      {
        const char* text = "Can't re-meshed a mesh with coincident free edges";
        SMESH_BadInputElements* error =
          new SMESH_BadInputElements( myMeshDS, COMPERR_BAD_INPUT_MESH, text );
        for ( size_t i = 0; i < freeBords._borders.size(); ++i )
          error->myBadElements.insert( error->myBadElements.end(),
                                       freeBords._borders[i].begin(),
                                       freeBords._borders[i].end() );
        theMesh.GetSubMesh( theMesh.GetShapeToMesh() )->GetComputeError().reset( error );
        throw SALOME_Exception( text );
      }
    }

    // fill holes
    myHole.resize( holes.size() );
    myInHolePos.resize( holes.size() );
    std::vector<const SMDS_MeshElement*> newFaces;
    for ( size_t i = 0; i < holes.size(); ++i )
    {
      newFaces.clear();
      SMESH_MeshAlgos::FillHole( holes[i], *myMeshDS, newFaces );

      // keep data to be able to remove hole filling faces after remeshing
      if ( !newFaces.empty() )
      {
        myHole[i].resize( holes[i].size() );
        for ( size_t iP = 0; iP < holes[i].size(); ++iP )
          myHole[i][iP] = SMESH_NodeXYZ( holes[i][iP] );

        myInHolePos[i] = ( SMESH_NodeXYZ( newFaces[0]->GetNode(0)) +
                           SMESH_NodeXYZ( newFaces[0]->GetNode(1)) +
                           SMESH_NodeXYZ( newFaces[0]->GetNode(2)) ) / 3.;
        myCapElems.insert( myCapElems.end(), newFaces.begin(), newFaces.end() );
        // unmark to be able to remove them if meshing is canceled
        // for ( size_t iF = 0; iF < newFaces.size(); ++iF )
        //   newFaces[iF]->setIsMarked( false );
      }
    }
    // fix orientation
    if ( !isGoodOri )
    {
      SMDS_ElemIteratorPtr fIt = myMeshDS->elementsIterator( SMDSAbs_Face );
      while ( fIt->more() )
      {
        const SMDS_MeshElement* f = fIt->next();
        gp_XYZ normal;
        if ( SMESH_MeshAlgos::FaceNormal( f, normal ))
        {
          TIDSortedElemSet allFaces;
          editor.Reorient2D( allFaces, normal, f );
          break;
        }
      }
    }
  }
  //================================================================================
  /*!
   * \brief Add hole borders to be kept in a new mesh
   */
  //================================================================================

  void HoleFiller::AddHoleBordersAndEdges( Ng_STL_Geometry * ngStlGeo, bool toAddEdges )
  {
    nglib::readedges.SetSize(0);

    for ( size_t i = 0; i < myHole.size(); ++i )
      for ( size_t iP = 1; iP < myHole[i].size(); ++iP )
      {
        Ng_STL_AddEdge( ngStlGeo,
                        myHole[i][iP-1].ChangeData(),
                        myHole[i][iP-0].ChangeData() );
      }

    if ( toAddEdges )
    {
      std::vector<const SMDS_MeshNode *>    nodes(2);
      std::vector<const SMDS_MeshElement *> faces(2);
      SMDS_EdgeIteratorPtr eIt = myMeshDS->edgesIterator();
      while ( eIt->more() )
      {
        const SMDS_MeshElement* edge = eIt->next();
        nodes[0] = edge->GetNode(0);
        nodes[1] = edge->GetNode(1);
        // check that an edge is a face border
        if ( myMeshDS->GetElementsByNodes( nodes, faces, SMDSAbs_Face ))
        {
          Ng_STL_AddEdge( ngStlGeo,
                          SMESH_NodeXYZ( nodes[0] ).ChangeData(),
                          SMESH_NodeXYZ( nodes[1] ).ChangeData() );
        }
      }
    }
    return;
  }
  //================================================================================
  /*!
   * \brief Remove triangles filling the holes
   */
  //================================================================================

  HoleFiller::~HoleFiller()
  {
    if ( myMeshDS->NbNodes() < 3 )
      return;

    if ( !myCapElems.empty() ) // old mesh not removed; simply remove myCapElems
    {
      for ( size_t i = 0; i < myCapElems.size(); ++i )
        myMeshDS->RemoveFreeElement( myCapElems[i], /*sm=*/0 );
      return;
    }

    bool hasOrphanNodes = true;

    const double tol = Max( 1e-3 * netgen::mparam.minh, Precision::Confusion() );

    for ( size_t i = 0; i < myHole.size(); ++i )
    {
      std::vector< gp_XYZ >& borderPnt = myHole[i];
      const gp_XYZ&          inHolePos = myInHolePos[i];
      if ( borderPnt.empty() ) continue;
      borderPnt.pop_back(); // first point repeated at end

      // mark all nodes located on the hole border

      // new nodeSearcher for each hole, otherwise it contains removed nodes for i > 0
      SMESHUtils::Deleter< SMESH_NodeSearcher > nodeSearcher;
      if ( hasOrphanNodes )
      {
        std::vector< const SMDS_MeshNode* > sharedNodes;
        sharedNodes.reserve( myMeshDS->NbNodes() );
        SMDS_NodeIteratorPtr nIt = myMeshDS->nodesIterator();
        while ( nIt->more() )
        {
          const SMDS_MeshNode* n = nIt->next();
          if ( n->NbInverseElements() )
            sharedNodes.push_back( n );
        }
        hasOrphanNodes = ((int) sharedNodes.size() < myMeshDS->NbNodes() );
        SMDS_ElemIteratorPtr elemIt( new SMDS_NodeVectorElemIterator( sharedNodes.begin(),
                                                                      sharedNodes.end() ));
        nodeSearcher._obj = SMESH_MeshAlgos::GetNodeSearcher( elemIt );
      }
      else
      {
        nodeSearcher._obj = SMESH_MeshAlgos::GetNodeSearcher( *myMeshDS );
      }

      std::vector< const SMDS_MeshElement* > edgesToRemove;
      edgesToRemove.reserve( borderPnt.size() );

      // look for a border point coincident with a node
      size_t iP = 0;
      SMESH_NodeXYZ bordNode1;
      for ( ; iP < borderPnt.size(); ++iP )
      {
        bordNode1 = nodeSearcher->FindClosestTo( borderPnt[iP] );
        if (( bordNode1 - borderPnt[iP] ).SquareModulus() < tol * tol )
          break;
      }
      ++iP;
      bordNode1._node->setIsMarked( true );

      // find the rest nodes located on the hole border
      boost::container::flat_set< const SMDS_MeshNode* > checkedNodes;
      gp_XYZ p1 = bordNode1;
      for ( size_t j = 0; j < borderPnt.size()+1; ++j,  iP = ( iP+1 ) % borderPnt.size() )
      {
        // among nodes surrounding bordNode1 find one most close to vec12
        gp_XYZ vec12 = borderPnt[iP] - p1;
        bool pntReached = false; // last found node is at iP
        while ( !pntReached )
        {
          const SMDS_MeshNode* bordNode = bordNode1._node;
          SMDS_ElemIteratorPtr fIt = bordNode->GetInverseElementIterator( SMDSAbs_Face );
          double minArea = 1e100;
          checkedNodes.clear();
          checkedNodes.insert( bordNode );
          while ( fIt->more() )
          {
            const SMDS_MeshElement* f = fIt->next();
            for ( int iN = 0, nbN = f->NbNodes(); iN < nbN; ++iN )
            {
              const SMDS_MeshNode* n = f->GetNode( iN );
              if ( !checkedNodes.insert( n ).second )
                continue;
              SMESH_NodeXYZ pn = n;
              gp_XYZ vecPN = pn - bordNode1;
              if ( vecPN * vec12 <= 0 )
                continue;
              gp_XYZ vec1N = pn - p1;
              double     a = vec12.CrossSquareMagnitude( vec1N );
              if ( a < minArea )
              {
                bordNode = n;
                minArea = a;
              }
            }
            if ( minArea < std::numeric_limits<double>::min() )
              break;
          }
          if ( bordNode == bordNode1._node )
            return; // bug in the loop above

          SMESH_NodeXYZ bordNode2 = bordNode;
          gp_XYZ            vec1N = bordNode2 - p1;
          double u = ( vec12 * vec1N ) / vec12.SquareModulus(); // param [0,1] of bordNode on vec12
          if ( u < 1 + tol )
          {
            bordNode->setIsMarked( true );
            //cout << bordNode->GetID() << " ";

            if ( const SMDS_MeshElement* edge = myMeshDS->FindEdge( bordNode1._node, bordNode ))
              edgesToRemove.push_back( edge );
            else
              edgesToRemove.push_back( myMeshDS->AddEdge( bordNode1._node, bordNode ));
            edgesToRemove.back()->setIsMarked( true );

            if ( minArea > std::numeric_limits<double>::min() &&
                 minArea / vec12.SquareModulus() > tol * tol )
            {
              // node is far from the border, move it
              gp_XYZ p = p1 + u * vec12;
              myMeshDS->MoveNode( bordNode, p.X(), p.Y(), p.Z() );
            }
            bordNode1 = bordNode2;
          }
          //else -- there must be another border point between bordNode1 and bordNode
          pntReached = ( u > 1 - tol );
        }
        p1 = borderPnt[iP];

      }
      //cout << endl << endl;

      // remove all faces starting from inHolePos

      // get a starting face
      std::vector< const SMDS_MeshNode* >     nodesToRemove;
      std::vector< const SMDS_MeshElement* >  facesToRemove;
      const SMDS_MeshNode* inHoleNode = nodeSearcher->FindClosestTo( inHolePos );
      if ( inHoleNode && ! inHoleNode->isMarked() )
      {
        SMDS_ElemIteratorPtr fIt = inHoleNode->GetInverseElementIterator( SMDSAbs_Face );
        while ( fIt->more() )
          facesToRemove.push_back( fIt->next() );
      }
      else
      {
        SMESHUtils::Deleter< SMESH_ElementSearcher > faceSearcher
          ( SMESH_MeshAlgos::GetElementSearcher( *myMeshDS ));
        if ( const SMDS_MeshElement* f = faceSearcher->FindClosestTo( inHolePos, SMDSAbs_Face ))
          facesToRemove.push_back( f );
        else
          continue;
      }
      for ( size_t iF = 0; iF < facesToRemove.size(); ++iF )
        facesToRemove[iF]->setIsMarked( true );

      // remove faces and nodes
      TIDSortedElemSet elemSet, avoidSet;
      const SMDS_MeshElement* e;
      while ( !facesToRemove.empty() )
      {
        const SMDS_MeshElement* inHoleFace = facesToRemove.back();
        facesToRemove.pop_back();

        // add adjacent faces into facesToRemove
        for ( int iN = 0, nbN = inHoleFace->NbNodes(); iN < nbN; ++iN )
        {
          const SMDS_MeshNode* n1 = inHoleFace->GetNode( iN );
          if ( !n1->isMarked() )
          {
            SMDS_ElemIteratorPtr eIt = n1->GetInverseElementIterator();
            while ( eIt->more() )
            {
              e = eIt->next();
              if ( e->GetType() == SMDSAbs_Face )
              {
                if ( !e->isMarked() )
                  facesToRemove.push_back( e );
                e->setIsMarked( true );
              }
              else if ( e->GetType() == SMDSAbs_Edge )
              {
                myMeshDS->RemoveFreeElement( e, 0, /*fromGroups=*/false );
              }
            }
            if ( n1->NbInverseElements() == 1 )
              nodesToRemove.push_back( n1 );
          }
          else
          {
            const SMDS_MeshNode* n2 = inHoleFace->GetNodeWrap( iN+1 );
            if (( n2->isMarked() ) &&
                ( !(e = myMeshDS->FindEdge( n1, n2 )) || !e->isMarked() )) // n1-n2 not hole border
            {
              if ( e ) // remove edge
                myMeshDS->RemoveFreeElement( e, 0, /*fromGroups=*/false );
              avoidSet.clear();
              avoidSet.insert( inHoleFace );
              if (( e = SMESH_MeshAlgos::FindFaceInSet( n1, n2, elemSet, avoidSet )))
              {
                if ( !e->isMarked() )
                  facesToRemove.push_back( e );
                e->setIsMarked( true );
              }
            }
          }
        }
        myMeshDS->RemoveFreeElement( inHoleFace, 0, /*fromGroups=*/false );

        for ( size_t iN = 0; iN < nodesToRemove.size(); ++iN )
          myMeshDS->RemoveFreeNode( nodesToRemove[iN], 0, /*fromGroups=*/false );
        nodesToRemove.clear();
      }

      // remove edges from the hole border
      // for ( size_t iE = 0; iE < edgesToRemove.size(); ++iE )
      //   myMeshDS->RemoveFreeElement( edgesToRemove[iE], 0, /*fromGroups=*/false );

    } // loop on holes

    return;
  } // ~HoleFiller()


  //================================================================================
  /*!
   * \brief Fix nodes of a group
   */
  //================================================================================

  void fixNodes( SMESHDS_GroupBase* group, netgen::STLGeometry* stlGeom )
  {
    SMESH_MeshAlgos::MarkElemNodes( group->GetElements(), false ); // un-mark nodes

    for ( SMDS_ElemIteratorPtr eIt = group->GetElements(); eIt->more(); )
    {
      const SMDS_MeshElement* e = eIt->next();
      for ( SMDS_NodeIteratorPtr nIt = e->nodeIterator(); nIt->more(); )
      {
        const SMDS_MeshNode* n = nIt->next();
        if ( n->isMarked() )
          continue;
        n->setIsMarked( true );

        SMESH_NodeXYZ p( n );
        int id = stlGeom->GetPointNum( netgen::Point<3>( p.X(),p.Y(),p.Z() ));
        if ( id > 0 )
          stlGeom->SetLineEndPoint( id );
      }
    }
  }

} // namespace

//=============================================================================
/*!
 * Constructor
 */
//=============================================================================

NETGENPlugin_Remesher_2D::NETGENPlugin_Remesher_2D(int hypId, SMESH_Gen* gen)
  : SMESH_2D_Algo(hypId, gen)
{
  _name = "NETGEN_Remesher_2D";
  _shapeType = (1 << TopAbs_FACE); // 1 bit /shape type
  _compatibleHypothesis.push_back("NETGEN_RemesherParameters_2D");
  _requireShape = false;

  _hypothesis = 0;
}

//=============================================================================
/*!
 * Check assigned hypotheses
 */
//=============================================================================

bool NETGENPlugin_Remesher_2D::CheckHypothesis (SMESH_Mesh&         theMesh,
                                                const TopoDS_Shape& theShape,
                                                Hypothesis_Status&  theStatus)
{
  _hypothesis = 0;

  // can work with no hypothesis
  theStatus = SMESH_Hypothesis::HYP_OK;

  const list<const SMESHDS_Hypothesis*>& hyps =
    GetUsedHypothesis( theMesh, theShape, /*skipAux=*/true );

  switch ( hyps.size() ) {
  case 0:
    break;
  case 1:
    _hypothesis = hyps.front();
    break;
  default:
    theStatus = SMESH_Hypothesis::HYP_INCOMPATIBLE;
  }

  return theStatus == SMESH_Hypothesis::HYP_OK;
}

//=============================================================================
/*!
 * Compute mesh on an input mesh
 */
//=============================================================================

bool NETGENPlugin_Remesher_2D::Compute(SMESH_Mesh&         theMesh,
                                       SMESH_MesherHelper* theHelper)
{
  if ( theMesh.NbFaces() == 0 )
    return !error( COMPERR_WARNING, "No faces in input mesh");

  NETGENPlugin_Mesher mesher( &theMesh, theMesh.GetShapeToMesh(), /*isVol=*/false);
  NETGENPlugin_NetgenLibWrapper ngLib;
  netgen::Mesh *        ngMesh = (netgen::Mesh*) ngLib._ngMesh;
  Ng_STL_Geometry *   ngStlGeo = Ng_STL_NewGeometry();
  netgen::STLTopology* stlTopo = (netgen::STLTopology*) ngStlGeo;
  netgen::multithread.terminate = 0;

  const NETGENPlugin_RemesherHypothesis_2D* hyp =
    dynamic_cast<const NETGENPlugin_RemesherHypothesis_2D*>( _hypothesis );
  mesher.SetParameters( hyp );// for holeFiller

  SMESHDS_Mesh* meshDS = theMesh.GetMeshDS();
  HoleFiller holeFiller( theMesh );
  //theHelper->SetIsQuadratic( theMesh.NbFaces( ORDER_QUADRATIC ));

  // fill ngStlGeo with triangles
  SMDS_ElemIteratorPtr fIt = meshDS->elementsIterator( SMDSAbs_Face );
  while ( fIt->more() )
  {
    const SMDS_MeshElement* f = fIt->next();
    SMESH_NodeXYZ n1 = f->GetNode( 0 );
    SMESH_NodeXYZ n2 = f->GetNode( 1 );
    SMESH_NodeXYZ n3 = f->GetNode( 2 );
    Ng_STL_AddTriangle( ngStlGeo,
                        n1.ChangeData(),
                        n2.ChangeData(),
                        n3.ChangeData() );
    if ( f->NbNodes() > 3 )
    {
      n2.Set( f->GetNode( 3 ));
      Ng_STL_AddTriangle( ngStlGeo,
                          n1.ChangeData(),
                          n3.ChangeData(),
                          n2.ChangeData());
    }
  }
  // add edges
  bool toAddExistingEdges = ( hyp && hyp->GetKeepExistingEdges() );
  holeFiller.AddHoleBordersAndEdges( ngStlGeo, toAddExistingEdges );

  // init stl DS
  //netgen::stldoctor.geom_tol_fact = 1e-12; // pointtol=boundingbox.Diam()*stldoctor.geom_tol_fact
  Ng_Result ng_res = Ng_STL_InitSTLGeometry( ngStlGeo );
  if ( ng_res != NG_OK )
  {
#ifdef _DEBUG_
    holeFiller.KeepHole();
#endif
    std::string txt = "Error Initialising the STL Geometry";
    if ( !stlTopo->GetStatusText().empty() )
      txt += ". " + stlTopo->GetStatusText();
    return error( COMPERR_BAD_INPUT_MESH, txt );
  }

  // set parameters
  Ng_Meshing_Parameters ngParams;
  if ( hyp )
  {
    ngParams.maxh              = hyp->GetMaxSize();
    ngParams.minh              = hyp->GetMinSize();
    ngParams.meshsize_filename = (char*) hyp->GetMeshSizeFile().c_str();
    ngParams.quad_dominated    = hyp->GetQuadAllowed();

    netgen::stlparam.yangle                  = hyp->GetRidgeAngle();
    netgen::stlparam.edgecornerangle         = hyp->GetEdgeCornerAngle();
    netgen::stlparam.chartangle              = hyp->GetChartAngle();
    netgen::stlparam.outerchartangle         = hyp->GetOuterChartAngle();
    netgen::stlparam.resthchartdistfac       = hyp->GetRestHChartDistFactor();
    netgen::stlparam.resthchartdistenable    = hyp->GetRestHChartDistEnable();
    netgen::stlparam.resthlinelengthfac      = hyp->GetRestHLineLengthFactor();
    netgen::stlparam.resthlinelengthenable   = hyp->GetRestHLineLengthEnable();
#ifndef NETGEN_V6
    netgen::stlparam.resthcloseedgefac       = hyp->GetRestHCloseEdgeFactor();
    netgen::stlparam.resthcloseedgeenable    = hyp->GetRestHCloseEdgeEnable();
#endif
    netgen::stlparam.resthsurfcurvfac        = hyp->GetRestHSurfCurvFactor();
    netgen::stlparam.resthsurfcurvenable     = hyp->GetRestHSurfCurvEnable();
    netgen::stlparam.resthedgeanglefac       = hyp->GetRestHEdgeAngleFactor();
    netgen::stlparam.resthedgeangleenable    = hyp->GetRestHEdgeAngleEnable();
    netgen::stlparam.resthsurfmeshcurvfac    = hyp->GetRestHSurfMeshCurvFactor();
    netgen::stlparam.resthsurfmeshcurvenable = hyp->GetRestHSurfMeshCurvEnable();

    mesher.SetParameters( hyp );
  }
  else
  {
    double diagSize = Dist( stlTopo->GetBoundingBox().PMin(), stlTopo->GetBoundingBox().PMax());
    netgen::mparam.maxh = diagSize / GetGen()->GetBoundaryBoxSegmentation();
    netgen::mparam.minh = netgen::mparam.maxh;
  }

  // save netgen::mparam as Ng_STL_MakeEdges() modify it by Ng_Meshing_Parameters
  netgen::MeshingParameters savedParams = netgen::mparam;

  if ( SMESH_Group* fixedEdges = ( hyp ? hyp->GetFixedEdgeGroup( theMesh ) : 0 ))
  {
    netgen::STLGeometry* stlGeom = (netgen::STLGeometry*)ngStlGeo;

    // the following code is taken from STLMeshing() method
#ifdef NETGEN_V6
    stlGeom->Clear();
    stlGeom->BuildEdges( netgen::stlparam );
    stlGeom->MakeAtlas( *ngMesh, netgen::mparam, netgen::stlparam );
    stlGeom->CalcFaceNums();
    stlGeom->AddFaceEdges();
    fixNodes( fixedEdges->GetGroupDS(), stlGeom );
    stlGeom->LinkEdges( netgen::stlparam );
#else
    stlGeom->Clear();
    stlGeom->BuildEdges();
    stlGeom->MakeAtlas( *ngMesh );
    stlGeom->CalcFaceNums();
    stlGeom->AddFaceEdges();
    fixNodes( fixedEdges->GetGroupDS(), stlGeom );
    stlGeom->LinkEdges();
#endif
    ngMesh->ClearFaceDescriptors();
    for (int i = 1; i <= stlGeom->GetNOFaces(); i++)
      ngMesh->AddFaceDescriptor (netgen::FaceDescriptor (i, 1, 0, 0));

    stlGeom->edgesfound = 1;
  }
  else
  {
    Ng_STL_MakeEdges( ngStlGeo, ngLib.ngMesh(), &ngParams );
  }

  netgen::mparam = savedParams;

  double h = netgen::mparam.maxh;
  ngMesh->SetGlobalH( h );
  ngMesh->SetMinimalH( netgen::mparam.minh );
  ngMesh->SetLocalH( stlTopo->GetBoundingBox().PMin() - netgen::Vec3d(h, h, h),
                     stlTopo->GetBoundingBox().PMax() + netgen::Vec3d(h, h, h),
                     netgen::mparam.grading );
  ngMesh->LoadLocalMeshSize( ngParams.meshsize_filename );

  netgen::OCCGeometry occgeo;
  mesher.SetLocalSize( occgeo, *ngMesh );

  // meshing
  try
  {
    ng_res = Ng_STL_GenerateSurfaceMesh( ngStlGeo, ngLib.ngMesh(), &ngParams );
  }
  catch (netgen::NgException & ex)
  {
    if ( netgen::multithread.terminate )
      if ( !hyp || !hyp->GetLoadMeshOnCancel() )
        return false;
  }
  if ( ng_res != NG_OK )
    return error( "Error in Surface Meshing" );

  int nbN = ngMesh->GetNP();
  int nbE = ngMesh->GetNSeg();
  int nbF = ngMesh->GetNSE();
  if ( nbF == 0 )
    return error( "Error in Surface Meshing" );

  // remove existing mesh
  holeFiller.ClearCapElements();
  SMDS_ElemIteratorPtr eIt = meshDS->elementsIterator();
  while ( eIt->more() )
    meshDS->RemoveFreeElement( eIt->next(), /*sm=*/0 );
  SMDS_NodeIteratorPtr nIt = meshDS->nodesIterator();
  while ( nIt->more() )
    meshDS->RemoveFreeNode( nIt->next(), /*sm=*/0 );

  // retrieve new mesh

  // add nodes
  std::vector< const SMDS_MeshNode* > newNodes( nbN+1 );
  for ( int i = 1; i <= nbN; ++i )
  {
    const netgen::MeshPoint& p = ngMesh->Point(i);
    newNodes[i] = meshDS->AddNode( p(0),p(1),p(2) );
  }

  // add edges
  std::vector<const SMDS_MeshNode*> nodes(4);
  for ( int i = 1; i <= nbE; ++i )
  {
    const netgen::Segment& seg = ngMesh->LineSegment(i);
    nodes.clear();
    for ( int j = 0; j < 2; ++j )
    {
      size_t pind = seg.pnums[j];
      if ( pind > 0 && pind < newNodes.size() )
        nodes.push_back( newNodes[ pind ]);
      else
        break;
    }
    if ( nodes.size() == 2 && !meshDS->FindEdge( nodes[0], nodes[1] ))
      meshDS->AddEdge( nodes[0], nodes[1] );
  }

  // find existing groups
  const char* theNamePrefix = "Surface_";
  const int   theNamePrefixLen = strlen( theNamePrefix );
  std::vector< SMESHDS_Group* > groups;
  if ( hyp && hyp->GetMakeGroupsOfSurfaces() )
  {
    SMESH_Mesh::GroupIteratorPtr grIt = theMesh.GetGroups();
    while ( grIt->more() )
    {
      SMESH_Group* group = grIt->next();
      SMESHDS_Group* groupDS;
      if (( group->GetGroupDS()->GetType() == SMDSAbs_Face ) &&
          ( strncmp( group->GetName(), theNamePrefix, theNamePrefixLen ) == 0 ) &&
          ( groupDS = dynamic_cast<SMESHDS_Group*>( group->GetGroupDS() )))
        groups.push_back( groupDS );
    }
  }

  // add faces
  for ( int i = 1; i <= nbF; ++i )
  {
    const netgen::Element2d& elem = ngMesh->SurfaceElement(i);
    nodes.clear();
    for ( int j = 1; j <= elem.GetNP(); ++j )
    {
      size_t pind = elem.PNum(j);
      if ( pind > 0 && pind < newNodes.size() )
        nodes.push_back( newNodes[ pind ]);
      else
        break;
    }
    const SMDS_MeshElement* newFace = 0;
    switch( nodes.size() )
    {
    case 3: newFace = meshDS->AddFace( nodes[0], nodes[1], nodes[2] ); break;
    case 4: newFace = meshDS->AddFace( nodes[0], nodes[1], nodes[2], nodes[3] ); break;
    }

    // add newFace to a group
    if ( newFace && hyp && hyp->GetMakeGroupsOfSurfaces() )
    {
      if ((size_t) elem.GetIndex()-1 >= groups.size() )
        groups.resize( elem.GetIndex(), 0 );

      SMESHDS_Group* & group = groups[  elem.GetIndex()-1 ];
      if ( !group )
      {
        SMESH_Group* gr = theMesh.AddGroup( SMDSAbs_Face, "");
        group = static_cast<SMESHDS_Group*>( gr->GetGroupDS() );
      }
      group->SMDSGroup().Add( newFace );
    }
  }

  // update groups
  int groupIndex = 1;
  for ( size_t i = 0; i < groups.size(); ++i )
  {
    if ( !groups[i] )
      continue;
    if ( groups[i]->IsEmpty() )
    {
      theMesh.RemoveGroup( groups[i]->GetID() );
    }
    else if ( SMESH_Group* g = theMesh.GetGroup( groups[i]->GetID() ))
    {
      g->SetName( SMESH_Comment( theNamePrefix ) << groupIndex++ );
    }
  }

  // as we don't assign the new triangles to a shape (the pseudo-shape),
  // to avoid their removal at hypothesis modification,
  // we mark the shape as always computed to avoid the error messages
  // that no elements assigned to the shape
  theMesh.GetSubMesh( theHelper->GetSubShape() )->SetIsAlwaysComputed( true );

  return true;
}

//=============================================================================
/*!
 * Do not compute mesh on geometry
 */
//=============================================================================

bool NETGENPlugin_Remesher_2D::Compute(SMESH_Mesh&         /*theMesh*/,
                                       const TopoDS_Shape& /*theShape*/)
{
  return false;
}

//=============================================================================
/*!
 * Terminate Compute()
 */
//=============================================================================

void NETGENPlugin_Remesher_2D::CancelCompute()
{
  SMESH_Algo::CancelCompute();
  netgen::multithread.terminate = 1;
}

//================================================================================
/*!
 * \brief Return progress of Compute() [0.,1]
 */
//================================================================================

double NETGENPlugin_Remesher_2D::GetProgress() const
{
  return netgen::multithread.percent / 100.;
}

//=============================================================================
/*!
 *
 */
//=============================================================================

bool NETGENPlugin_Remesher_2D::Evaluate(SMESH_Mesh&         /*aMesh*/,
                                        const TopoDS_Shape& /*aShape*/,
                                        MapShapeNbElems& /*aResMap*/)
{
  return false;
}
