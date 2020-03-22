/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

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
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */

/*! \file MsqMeshEntity.hpp

\author Darryl Melander
\author Thomas Leurent
\author Michael Brewer

*/

#ifndef MSQMESHENTITY_HPP
#define MSQMESHENTITY_HPP

#include "Mesquite.hpp"
#include "TopologyInfo.hpp"
#include "Vector3D.hpp"
#include "NodeSet.hpp"
#include "Sample.hpp"

#include <vector>
#include <cstring>
#include <cassert>


namespace MESQUITE_NS
{
  class PatchData;
  class MsqVertex;

  /*!
      \class MsqMeshEntity
      \brief MsqMeshEntity is the Mesquite object that stores information about
      the elements in the mesh.

      
  */
  class MESQUITE_EXPORT MsqMeshEntity
  {
  public:
    
    MsqMeshEntity()
      : mType(MIXED), vertexIndices(0), numVertexIndices(0)
      {}

      //! Returns element type
    inline EntityTopology get_element_type() const
      { return mType; }
    
      //! Returns the number of vertices in this element,
      //! based on its element type.
    inline std::size_t vertex_count() const;
      //! Return number of nodes in element (number of corner
      //! vertices + number of higher-order nodes).
    inline std::size_t node_count() const 
      { return numVertexIndices; }
      //! Returns number of target matrices for this element type
    inline std::size_t corner_count() const
      { return mType == PYRAMID ? 4 : vertex_count(); }
    
      //! gets the vertices of the mesh entity
    void get_vertex_indices(std::vector<std::size_t> &vertex_list) const;
    void append_vertex_indices(std::vector<std::size_t> &vertex_list) const;
    size_t get_local_matrix_map_about_vertex(PatchData &pd,
                                             MsqVertex* vert,
                                             size_t local_map_size,
                                             int* local_map,
                                             MsqError &err) const;
      //! gets the vertices of the mesh entity
    void get_node_indices(std::vector<std::size_t> &vertex_list) const;
    void append_node_indices(std::vector<std::size_t> &vertex_list) const;
    //! Very efficient retrieval of vertices indexes 
    //! (corresponding to the PatchData vertex array).
    inline const std::size_t *get_vertex_index_array() const;
    inline std::size_t *get_vertex_index_array();
    
      //! Sets element data
    void set_element_type(EntityTopology type)
      { mType = type; }

      //! Set connectivity data (vertex array) for element.
      //! MsqMeshEntity keeps the pointer to the passed array, it
      //! does not copy it.  The caller is still responsible for
      //! releasing the memory for the passed index array after
      //! the MsqMeshEntity is destroyed.  The intention is that
      //! this is a pointer to a portion of a larger connectivity array
      //! managed by the owning PatchData.
    void set_connectivity( std::size_t *indices, size_t num_vertices);

    std::size_t get_vertex_index(std::size_t vertex_in_element) const;
    
    //! Returns the centroid of the element.
    void get_centroid(Vector3D& centroid, const PatchData &pd, MsqError &err) const;
    
      //!Fills a vector<size_t> with vertices connected to the given
      //!vertex through the edges of this MsqMeshEntity.
    void get_connected_vertices(std::size_t vertex_index,
                                std::vector<std::size_t> &vert_indices,
                                MsqError &err);
    
      //!Computes the area of the element.
      //!The returned value is always non-negative.
    double compute_unsigned_area(PatchData &pd, MsqError &err );
    
      //!Computes the signed area of the element.
    double compute_signed_area(PatchData &pd, MsqError &err );
    
      //! Check sign of the determinant mapping fuction Jacobian 
      //! at representative sample points in the element.
      //!\param inverted_count  Number of sampling locations within the
      //!                       element for which the determinant of the
      //!                       mapping function Jacobian was negative.
      //!\param tested_count    The number of sampling locations in the
      //!                       element for which the Jacobian was tested.
    void check_element_orientation( PatchData &pd, 
                                    int& inverted_count,
                                    int& tested_count,
                                    MsqError &err );
    void check_element_orientation_corners( PatchData &pd, 
                                            int& inverted_count,
                                            int& tested_count,
                                            MsqError &err );
    
    
      //! Uses a MeshDomain call-back function to compute the normal at the corner.
    //void compute_corner_normal( std::size_t corner_pt, 
    //                            Vector3D &normal,
    //                            PatchData &pd, 
    //                            MsqError &err);
    void compute_corner_normals( Vector3D normals[], PatchData& pd, MsqError& err );

      //! Get NodeSet indicating all nodes present for this element type.
    NodeSet all_nodes( MsqError& err ) const;

  private:

      //! Check for a negative Jacobian at the specified sample point
      //! of a volume element.
    bool inverted_jacobian_3d( PatchData& pd, NodeSet nodes, Sample sample, MsqError& err );

      //! Check for a negative Jacobian at the specified sample point
      //! of a surface element.
    bool inverted_jacobian_2d( PatchData& pd, NodeSet nodes, Sample sample, MsqError& err );
    
    EntityTopology mType;
    /** Pointer to connectivity array.
     *  NOTE: The memory occupied by this array is assumed to be
     *        owned/managed by the owning PatchData.  Do not try
     *        to delete/resize/etc. this array directly!
     */
    size_t* vertexIndices;  
    size_t numVertexIndices;
    
      // output operator for debugging.
    friend std::ostream& operator<<(std::ostream &stream, 
                                          const MsqMeshEntity &entity);
    
  };
  
  
    // Returns the number of vertices in this type
    // of element, or 0 if a variable number.
  inline size_t MsqMeshEntity::vertex_count() const
  { return mType == POLYGON || mType == POLYHEDRON ? 
           node_count() : TopologyInfo::corners( mType ); }

  inline void MsqMeshEntity::set_connectivity( std::size_t *indices,
                                               size_t num_vertices )
  {
    vertexIndices = indices;
    numVertexIndices = num_vertices;
  }
  
  inline const std::size_t *MsqMeshEntity::get_vertex_index_array() const
  { return vertexIndices; }
  
  inline std::size_t *MsqMeshEntity::get_vertex_index_array()
  { return vertexIndices; }
 
  inline std::size_t 
  MsqMeshEntity::get_vertex_index(std::size_t vertex_in_element) const
  {
      // Make sure we're in range
    assert(vertex_in_element < vertex_count());
      // Return the index
    return vertexIndices[vertex_in_element];
  }
} //namespace


#endif // MsqMeshEntity_hpp
