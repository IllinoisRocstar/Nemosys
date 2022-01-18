// Copyright (C) 2007-2020  CEA/DEN, EDF R&D, OPEN CASCADE
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

//  SMESH SMDS : implementation of Salome mesh data structure
//
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "SMDS_FaceOfNodes.hxx"

#include "SMDS_SetIterator.hxx"
#include "SMDS_MeshNode.hxx"
#include "SMDS_Mesh.hxx"

#include <utilities.h>

#include <boost/make_shared.hpp>

//=======================================================================
//function : NbEdges
//purpose  : 
//=======================================================================

int SMDS_FaceOfNodes::NbEdges() const
{
  return NbNodes();
}

int SMDS_FaceOfNodes::NbFaces() const
{
  return 1;
}

int SMDS_FaceOfNodes::NbNodes() const
{
  return myNbNodes;
}

int SMDS_FaceOfNodes::GetNodeIndex( const SMDS_MeshNode* node ) const
{
  for ( int i = 0; i < myNbNodes; ++i )
    if ( myNodes[i] == node )
      return i;
  return -1;
}

//=======================================================================
//function : Print
//purpose  : 
//=======================================================================

void SMDS_FaceOfNodes::Print(ostream & OS) const
{
  OS << "face <" << GetID() << " > : ";
  int i;
  for (i = 0; i < NbNodes() - 1; i++) OS << myNodes[i] << ",";
  OS << myNodes[i] << ") " << endl;
}

SMDS_ElemIteratorPtr SMDS_FaceOfNodes::nodesIterator() const
{
  return boost::make_shared< SMDS_NodeArrayElemIterator >( &myNodes[0], &myNodes[0] + NbNodes() );
}

SMDS_NodeIteratorPtr SMDS_FaceOfNodes::nodeIterator() const
{
  return boost::make_shared< SMDS_NodeArrayIterator >( &myNodes[0], &myNodes[0] + NbNodes() );
}

SMDS_FaceOfNodes::SMDS_FaceOfNodes(const SMDS_MeshNode* node1,
                                   const SMDS_MeshNode* node2,
                                   const SMDS_MeshNode* node3)
{
  myNbNodes = 3;
  myNodes[0]=node1;
  myNodes[1]=node2;
  myNodes[2]=node3;
  myNodes[3]=0;
}

SMDS_FaceOfNodes::SMDS_FaceOfNodes(const SMDS_MeshNode* node1,
                                   const SMDS_MeshNode* node2,
                                   const SMDS_MeshNode* node3,
                                   const SMDS_MeshNode* node4)
{
  myNbNodes = 4;
  myNodes[0]=node1;
  myNodes[1]=node2;
  myNodes[2]=node3;
  myNodes[3]=node4;       
}
bool SMDS_FaceOfNodes::ChangeNodes(const SMDS_MeshNode* nodes[],
                                   const int            nbNodes)
{
  myNbNodes = nbNodes;
  myNodes[0]=nodes[0];
  myNodes[1]=nodes[1];
  myNodes[2]=nodes[2];
  if (nbNodes == 4)
    myNodes[3]=nodes[3];
  else if (nbNodes != 3)
    return false;

  return true;
}

/*!
 * \brief Return node by its index
 * \param ind - node index
 * \retval const SMDS_MeshNode* - the node
 */
const SMDS_MeshNode* SMDS_FaceOfNodes::GetNode(const int ind) const
{
  return myNodes[ ind ];
}

SMDSAbs_EntityType SMDS_FaceOfNodes::GetEntityType() const
{
  return NbNodes() == 3 ? SMDSEntity_Triangle : SMDSEntity_Quadrangle;
}
SMDSAbs_GeometryType SMDS_FaceOfNodes::GetGeomType() const
{
  return NbNodes() == 3 ? SMDSGeom_TRIANGLE : SMDSGeom_QUADRANGLE;
}
