// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016-2020  EDF R&D
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
// See http://www.alneos.com/ or email : contact@alneos.fr
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef _GMSHPlugin_Hypothesis_2D_HXX_
#define _GMSHPlugin_Hypothesis_2D_HXX_

#include "GMSHPlugin_Defs.hxx"

#include "GMSHPlugin_Hypothesis.hxx"
#include "Utils_SALOME_Exception.hxx"

//  Parameters for work of GMSH.
// This class is just to give 2D dimension, actually
// it inherits all behaviour of the parent 

class GMSHPLUGIN_EXPORT  GMSHPlugin_Hypothesis_2D: public GMSHPlugin_Hypothesis
{
public:

  GMSHPlugin_Hypothesis_2D(int hypId, SMESH_Gen * gen);

};

#endif
