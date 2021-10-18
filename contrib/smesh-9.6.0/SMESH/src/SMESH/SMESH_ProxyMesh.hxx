// Copyright (C) 2007-2020  CEA/DEN, EDF R&D, OPEN CASCADE
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

// File      : SMESH_ProxyMesh.hxx
// Created   : Thu Dec  2 10:05:35 2010
// Author    : Edward AGAPOV (eap)

#ifndef __SMESH_ProxyMesh_HXX__
#define __SMESH_ProxyMesh_HXX__

#include "SMESH_SMESH.hxx"

#include "SMDS_ElementHolder.hxx"
#include "SMESHDS_SubMesh.hxx"
#include "SMESH_TypeDefs.hxx"

#include <TopoDS_Shape.hxx>
#include <NCollection_DataMap.hxx>

#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>

class SMDS_MeshNode;
class SMDS_MeshElement;
class SMESHDS_Mesh;
class SMESH_Mesh;

/*!
 * \brief Container of xD mesh elements substituting other ones in the
 *        input mesh of an (x+1)D algorithm
 */
class SMESH_EXPORT SMESH_ProxyMesh
{
public:

  typedef boost::shared_ptr<SMESH_ProxyMesh> Ptr;

  typedef std::map<const SMDS_MeshNode*, const SMDS_MeshNode*, TIDCompare > TN2NMap;

  //--------------------------------------------------------------------------------
  /*!
   * \brief Proxy sub-mesh
   */
  class SMESH_EXPORT SubMesh : public SMESHDS_SubMesh
  {
  public:

    const TN2NMap*       GetNodeNodeMap() const { return _n2n; }
    const SMDS_MeshNode* GetProxyNode( const SMDS_MeshNode* n ) const;
    const UVPtStructVec& GetUVPtStructVec() const { return _uvPtStructVec; }
    virtual void         AddElement(const SMDS_MeshElement * e);
    virtual int          NbElements() const;
    virtual int          NbNodes() const;
    virtual SMDS_ElemIteratorPtr GetElements() const;
    virtual SMDS_NodeIteratorPtr GetNodes() const;
    virtual void         Clear();
    virtual bool         Contains(const SMDS_MeshElement * ME) const;

    template< class ITERATOR >
      void ChangeElements( ITERATOR it, ITERATOR end )
    {
      // change SubMesh contents without deleting tmp elements
      // for which the caller is responsible
      _elements.assign( it, end );
    }
    SubMesh(const SMDS_Mesh* mesh, int index=0);
    virtual ~SubMesh() { Clear(); }

  protected:
    std::vector<const SMDS_MeshElement *> _elements;
    TN2NMap*                              _n2n;
    UVPtStructVec                         _uvPtStructVec; // for SubMesh of EDGE
    friend class SMESH_ProxyMesh;

  protected: // methods of SMDS_ElementHolder; remove elements before mesh compacting or clearing
    virtual SMDS_ElemIteratorPtr getElements() { Clear(); return GetElements(); }
    virtual void tmpClear() {}
    virtual void add( const SMDS_MeshElement* element ) {}
    virtual void compact() {}
    virtual void clear() { Clear(); }
  };
  //--------------------------------------------------------------------------------
  // Public interface

  SMESH_ProxyMesh();
  SMESH_ProxyMesh(std::vector<SMESH_ProxyMesh::Ptr>& components);
  SMESH_ProxyMesh(const SMESH_Mesh& mesh);
  virtual ~SMESH_ProxyMesh();

  // Returns a sub-mesh of a shape; it can be a proxy sub-mesh
  const SMESHDS_SubMesh* GetSubMesh(const TopoDS_Shape& shape) const;

  // Return a sub-mesh by a shape ID; it can be a proxy sub-mesh
  const SMESHDS_SubMesh* GetSubMesh(const int shapeID) const;

  // Return a proxy sub-mesh of a shape; it can be NULL
  const SubMesh*         GetProxySubMesh(const TopoDS_Shape& shape) const;

  // Return a proxy node of a node; the input node is returned if no proxy exists
  const SMDS_MeshNode*   GetProxyNode( const SMDS_MeshNode* node ) const;

  // Return number of proxy sub-meshes
  int                    NbProxySubMeshes() const;

  // Return iterator on all faces of the mesh taking into account substitutions.
  // To be used in case of mesh without shape
  SMDS_ElemIteratorPtr   GetFaces() const;

  // Return iterator on all faces on the face taking into account substitutions
  SMDS_ElemIteratorPtr   GetFaces(const TopoDS_Shape& face) const;

  // Return total nb of faces taking into account substitutions
  int                    NbFaces() const;

  bool                   IsTemporary(const SMDS_MeshElement* elem ) const;

  // Return iterator on inverse elements of a node that may be a proxy one
  SMDS_ElemIteratorPtr   GetInverseElementIterator(const SMDS_MeshNode* node,
                                                   SMDSAbs_ElementType  type) const;

  // Check if a FACE has prisms on its both sides
  static bool HasPrismsOnTwoSides( SMESHDS_SubMesh* faceSM );

  SMESH_Mesh*            GetMesh() const { return const_cast<SMESH_Mesh*>( _mesh ); }

  SMESHDS_Mesh*          GetMeshDS() const;

  //--------------------------------------------------------------------------------
  // Interface for descendants
 protected:

  void     setMesh(const SMESH_Mesh& mesh);

  int      shapeIndex(const TopoDS_Shape& shape) const;

  virtual SubMesh* newSubmesh(int index=0) const;

  // Return a proxy sub-mesh; zero index is for the case of mesh w/o shape
  SubMesh* findProxySubMesh(int shapeIndex=0) const;

  // Return a proxy sub-mesh; it is created if not yet exists
  SubMesh* getProxySubMesh(int shapeIndex);

  // Return a proxy sub-mesh; it is created if not yet exists
  SubMesh* getProxySubMesh(const TopoDS_Shape& shape=TopoDS_Shape());

  // move proxy sub-mesh from other proxy mesh to this, returns true if sub-mesh found
  bool     takeProxySubMesh( const TopoDS_Shape& shape, SMESH_ProxyMesh* proxyMesh );

  // move tmp elements residing the _mesh from other proxy mesh to this
  void     takeTmpElemsInMesh( SMESH_ProxyMesh* proxyMesh );

  // removes tmp element from the _mesh
  void     removeTmpElement( const SMDS_MeshElement* elem );

  // stores tmp element residing the _mesh
  void     storeTmpElement( const SMDS_MeshElement* elem );

  // store node-node correspondence
  void     setNode2Node(const SMDS_MeshNode* srcNode,
                        const SMDS_MeshNode* proxyNode,
                        const SubMesh*       subMesh);

  // types of elements needed to implement NbFaces() and GetFaces();
  // if _allowedTypes is empty, only elements from _subMeshes are returned,
  // else elements of _mesh filtered using allowedTypes are additionally returned
  std::vector< SMDSAbs_EntityType> _allowedTypes;

 private:

  const SMESH_Mesh*       _mesh;

  // proxy sub-meshes; index in vector == shapeIndex(shape)
  std::vector< SubMesh* > _subMeshes;

  // tmp elements residing the _mesh, to be deleted at destruction
  std::set< const SMDS_MeshElement* > _elemsInMesh;

  // inverse elements of proxy nodes
  typedef NCollection_DataMap< const SMDS_MeshNode*,
                               std::vector< const SMDS_MeshElement* >,
                               SMESH_Hasher                         > TNodeElemVecMap;
  TNodeElemVecMap  _inverseElements;

  // Complex submesh used to iterate over elements in other sub-meshes
  mutable SubMesh* _subContainer;
};

#endif
