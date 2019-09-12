// -*- C++ -*-
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

#ifndef _H_MAD_MADLIB
#define _H_MAD_MADLIB

#include "madlib_export.h"

// interface to the geometrical model
#include "ModelInterface.h"

// interface to the mesh
#include "MeshDataBaseInterface.h"
#include "MeshDataBaseParallelInterface.h"
#include "MeshDataBaseComm.h"
#include "MeshDataBaseCommPeriodic.h"
#include "CheckMesh.h"

// interface to mesh adaptation
#include "AdaptInterface.h"
#include "SizeField.h"
#include "CallbackDefinition.h"
#include "MobileObject.h"

//! Needed with the Gmsh geometrical model
void MADLIB_EXPORT MAdLibInitialize(int *argc, char **argv[]);
void MADLIB_EXPORT MAdLibFinalize();

#endif
