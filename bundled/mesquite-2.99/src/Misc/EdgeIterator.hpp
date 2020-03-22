/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file EdgeIterator.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_EDGE_ITERATOR_HPP
#define MSQ_EDGE_ITERATOR_HPP

#include "Mesquite.hpp"
#include "PatchData.hpp"
#include "MsqError.hpp"

namespace MESQUITE_NS {

/**\brief Iterate over all edges in a patch*/
class EdgeIterator
{
public: 
  EdgeIterator( PatchData* patch, MsqError& err );
  inline bool is_at_end() const;
  inline const Vector3D& start() const;
  inline const Vector3D& end() const;
  inline const Vector3D* mid() const;
  inline void step( MsqError& err );

  struct Edge {
    Edge( size_t vtx, size_t mid ) : otherVertex(vtx), midVertex(mid) {}
    Edge() {}
    size_t otherVertex;
    size_t midVertex;
  };
private:
  PatchData* patchPtr;
  size_t vertIdx;
  std::vector<Edge> adjList;
  std::vector<Edge>::iterator adjIter;
  void get_adjacent_vertices( MsqError& err );
};


inline bool operator<( const EdgeIterator::Edge& e1, const EdgeIterator::Edge& e2 )
  { return e1.otherVertex < e2.otherVertex || 
          (e1.otherVertex == e2.otherVertex && e1.midVertex < e2.midVertex); }

inline bool operator==( const EdgeIterator::Edge& e1, const EdgeIterator::Edge& e2 )
  { return e1.otherVertex == e2.otherVertex && e1.midVertex == e2.midVertex; }

bool EdgeIterator::is_at_end() const
{
  return vertIdx >= patchPtr->num_nodes();
}

const Vector3D& EdgeIterator::start() const
  { return patchPtr->vertex_by_index( vertIdx ); }

const Vector3D& EdgeIterator::end() const
  { return patchPtr->vertex_by_index( adjIter->otherVertex ); }

const Vector3D* EdgeIterator::mid() const
  { return adjIter->midVertex < patchPtr->num_nodes() ? 
           &patchPtr->vertex_by_index( adjIter->midVertex ) : 0; }

void EdgeIterator::step( MsqError& err )
{
  if (adjIter != adjList.end())
  {
    ++adjIter;
  }
  
  while (adjIter == adjList.end() && ++vertIdx < patchPtr->num_nodes())
  {
    get_adjacent_vertices( err );  MSQ_ERRRTN(err);
  }
}


} // namespace MESQUITE_NS

#endif
