/*
  Copyright (C) 2016-2017 The University of Edinburgh

  The file is part of libsupermesh
    
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation;
  version 2.1 of the License.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
 * The following code is derived from femtools/Halos_IO.cpp and
 * femtools/Tokenize.cpp in Fluidity git
 * revision 4e6c1d2b022df3a519cdec120fad28e60d1b08d9 (dated 2015-02-25)
 */

// Fluidity copyright information (note that AUTHORS mentioned in the following
// has been renamed to Fluidity_AUTHORS):

/*
  Copyright (C) 2006 Imperial College London and others.
  
  Please see the AUTHORS file in the main source directory for a full list
  of copyright holders.

  Prof. C Pain
  Applied Modelling and Computation Group
  Department of Earth Science and Engineering
  Imperial College London

  amcgsoftware@imperial.ac.uk
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation,
  version 2.1 of the License.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
  USA
*/

#include "tinyxml.h"

#include <cstddef>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "Read_Halos_C++.h"

using namespace std;
using namespace libsupermesh;

void libsupermesh::ReadHalos(const string &filename, int &process, int &nprocs,
  map<int, int> &npnodes, map<int, vector<vector<int> > > &send,
  map<int, vector<vector<int> > > &recv) { 
  npnodes.clear();
  send.clear();
  recv.clear();
  
  // Read the halo file
  TiXmlDocument doc(filename.c_str());
  if(!doc.LoadFile()) {
    cerr << doc.ErrorDesc() << endl;
    string error = string("Failed to read halo file '") + filename + string("'");
    libsupermesh_abort(error.c_str());
  }
  
  const char *buffer;
   
  // Extract the XML header
  TiXmlNode *header = doc.FirstChild();
  while(header and header->Type() != TiXmlNode::TINYXML_DECLARATION) {
    header = header->NextSibling();
  }
  if(!header) {
    string error = string("Invalid halo file '") + filename + string("': Missing XML declaration");
    libsupermesh_abort(error.c_str());
  }

  // Extract the root node
  TiXmlNode *rootNode = header->NextSiblingElement();
  if(!rootNode) {
    string error = string("Invalid halo file '") + filename + string("': Missing root element");
    libsupermesh_abort(error.c_str());
  }
  TiXmlElement *rootEle = rootNode->ToElement();
  
  // Extract process
  buffer = rootEle->Attribute("process");
  if(!buffer) {
    string error = string("Invalid halo file '") + filename + string("': Missing process attribute");
    libsupermesh_abort(error.c_str());
  }
  process = atoi(buffer);
  if(process < 0) {
    string error = string("Invalid halo file '") + filename + string("': Invalid process attribute");
    libsupermesh_abort(error.c_str());
  }
  
  // Extract nprocs
  buffer = rootEle->Attribute("nprocs");
  if(!buffer) {
    string error = string("Invalid halo file '") + filename + string("': Missing nprocs attribute");
    libsupermesh_abort(error.c_str());
  }
  nprocs = atoi(buffer);
  if(nprocs < 1) {
    string error = string("Invalid halo file '") + filename + string("': Invalid nprocs attribute");
    libsupermesh_abort(error.c_str());
  } else if(process >= nprocs) {
    string error = string("Invalid halo file '") + filename + string("': Invalid process / nprocs attributes");
    libsupermesh_abort(error.c_str());
  }
  
  // Extract halo data for each process for each level
  // Find the next halo element
  for(TiXmlNode *haloNode = rootEle->FirstChildElement("halo");haloNode;haloNode = haloNode->NextSiblingElement("halo")) {
    TiXmlElement *haloEle = haloNode->ToElement();
    
    // Extract the level
    buffer = haloEle->Attribute("level");
    if(!buffer) {
      string error = string("Invalid halo file '") + filename + string("': halo_data element missing level attribute");
      libsupermesh_abort(error.c_str());
    }
    // Check that data for this level has not already been extracted
    int level = atoi(buffer);
    if(send.count(level) > 0 or recv.count(level) > 0) {
      string error = string("Invalid halo file '") + filename + string("': Multiple halos defined for a single level");
      libsupermesh_abort(error.c_str());
    }
    send[level] = vector<vector<int> >(nprocs);
    recv[level] = vector<vector<int> >(nprocs);

    // Extract n_private_nodes
    buffer = haloEle->Attribute("n_private_nodes");
    if(!buffer) {      
      string error = string("Invalid halo file '") + filename + string("': halo_data element missing n_private_nodes attribute");
      libsupermesh_abort(error.c_str());
    }
    npnodes[level] = atoi(buffer);
    
    // Find the next halo_data element
    for(TiXmlNode *dataNode = haloEle->FirstChildElement("halo_data");dataNode;dataNode = dataNode->NextSiblingElement("halo_data")) {
      TiXmlElement *dataEle = dataNode->ToElement();
    
      // Extract the process
      buffer = dataEle->Attribute("process");
      if(!buffer) {
        string error = string("Invalid halo file '") + filename + string("': halo_data element missing process attribute");
        libsupermesh_abort(error.c_str());
      }
      int proc = atoi(buffer);
      if(proc < 0 or proc >= nprocs) {
        string error = string("Invalid halo file '") + filename + string("': Invalid halo_data element process / nprocs attributes");
        libsupermesh_abort(error.c_str());
      }
      
      // Check that data for this level and process has not already been extracted
      if(send[level][proc].size() > 0 or recv[level][proc].size() > 0) {        
        string error = string("Invalid halo file '") + filename + string("': Multiple halos defined for a single level and process");
        libsupermesh_abort(error.c_str());
      }
      
      // Extract the send data
      TiXmlNode *sendDataNode = dataEle->FirstChildElement("send");
      if(sendDataNode) {
        TiXmlNode *sendDataTextNode = sendDataNode->FirstChild();
        while(sendDataTextNode and sendDataTextNode->Type() != TiXmlNode::TINYXML_TEXT) {
          sendDataTextNode = sendDataTextNode->NextSibling();
        }
        if(sendDataTextNode) {
          vector<string> tokens;
          Tokenize(string(sendDataTextNode->Value()), tokens, " ");
          for(vector<string>::size_type i = 0;i < tokens.size();i++) {
            send[level][proc].push_back(atoi(tokens[i].c_str()));
          }
        }
        
        if(sendDataNode->NextSiblingElement("data")) {        
          string error = string("Invalid halo file '") + filename + string("': Multiple sets of sends defined for a single level and process");
          libsupermesh_abort(error.c_str());
        }
      }
      
      // Extract the receive data
      TiXmlNode *recvDataNode = dataEle->FirstChildElement("receive");
      if(recvDataNode) {
        TiXmlNode *recvDataTextNode = recvDataNode->FirstChild();
        while(recvDataTextNode and recvDataTextNode->Type() != TiXmlNode::TINYXML_TEXT) {
          recvDataTextNode = recvDataTextNode->NextSibling();
        }
        if(recvDataTextNode) {
          vector<string> tokens;
          Tokenize(string(recvDataTextNode->Value()), tokens, " ");
          for(vector<string>::size_type i = 0;i < tokens.size();i++) {
            recv[level][proc].push_back(atoi(tokens[i].c_str()));
          }
        }
        
        if(recvDataNode->NextSiblingElement("data")) {        
          string error = string("Invalid halo file '") + filename + string("': Multiple sets of receives defined for a single level and process");
          libsupermesh_abort(error.c_str());
        }
      }
    }
  }
}

void libsupermesh::Tokenize(const string &str, vector<string> &tokens, const string &delimiters) {
  tokens.clear();
  
  // Skip delimiter at beginning
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  
  // Find first delimiter
  string::size_type pos = str.find_first_of(delimiters, lastPos);
  
  while(lastPos != string::npos) {
    // Found a token, add it to the vector
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // Skip delimiters.  Note the "not_of".
    lastPos = str.find_first_not_of(delimiters, pos);

    // Find next delimiter
    pos = str.find_first_of(delimiters, lastPos);
  }
}
  
extern "C" {  
  void libsupermesh_read_halo(void **data, const char *filename, int *process,
    int *nprocs) {      
    HaloData *halo_data = new HaloData();
    *data = static_cast<void*>(halo_data);
  
    ReadHalos(string(filename),
      halo_data->process, halo_data->nprocs,
      halo_data->npnodes, halo_data->send, halo_data->recv);
    
    *process = halo_data->process;
    *nprocs = halo_data->nprocs;
  }
  
  void libsupermesh_halo_sizes(void **data, int level, int *nsends,
    int *nreceives) {
    HaloData *halo_data = static_cast<HaloData*>(*data);
    int nprocs = halo_data->nprocs;
    
    if(halo_data->npnodes.count(level) == 0 || 
       halo_data->send.count(level) == 0 ||
       halo_data->recv.count(level) == 0) {
       ostringstream error_buffer;
       error_buffer << "Data not found for halo level " << level << "";
       libsupermesh_abort(error_buffer.str().c_str());
     }
          
    for(int i = 0;i < nprocs;i++) {
      nsends[i] = halo_data->send[level][i].size();
      nreceives[i] = halo_data->recv[level][i].size();
    }
  }
  
  void libsupermesh_halo_data(void **data, int level, int *npnodes,
    int *send, int *recv) {
    HaloData *halo_data = static_cast<HaloData*>(*data);
    int nprocs = halo_data->nprocs;
    
    if(halo_data->npnodes.count(level) == 0 || 
       halo_data->send.count(level) == 0 ||
       halo_data->recv.count(level) == 0) {
       ostringstream error_buffer;
       error_buffer << "Data not found for halo level " << level << "";
       libsupermesh_abort(error_buffer.str().c_str());
     }
    
    *npnodes = halo_data->npnodes[level];
    
    ptrdiff_t index = 0;
    for(int i = 0;i < nprocs;i++) {
      std::vector<int>::size_type nsends = halo_data->send[level][i].size();
      memcpy(send + index, halo_data->send[level][i].data(), nsends * sizeof(int));
      index += nsends;
    }
    
    index = 0;
    for(int i = 0;i < nprocs;i++) {
      std::vector<int>::size_type nreceives = halo_data->recv[level][i].size();
      memcpy(recv + index, halo_data->recv[level][i].data(), nreceives * sizeof(int));
      index += nreceives;
    }
  }
    
  void libsupermesh_deallocate_halo(void **data) {
    delete (static_cast<HaloData*>(*data));
  }
}
