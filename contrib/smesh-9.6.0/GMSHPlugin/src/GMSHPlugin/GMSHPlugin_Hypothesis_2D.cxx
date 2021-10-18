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
#include "GMSHPlugin_Hypothesis_2D.hxx"
#include <utilities.h>

//#include <stdio.h>

using namespace std;

//=============================================================================
/*!
 *  
 */
//=============================================================================
GMSHPlugin_Hypothesis_2D::GMSHPlugin_Hypothesis_2D (int hypId, SMESH_Gen * gen)
  : GMSHPlugin_Hypothesis(hypId, gen)

{
  _name = "GMSH_Parameters_2D";
  _param_algo_dim = 2;
}

