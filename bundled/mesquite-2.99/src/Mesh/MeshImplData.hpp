/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */

#ifndef MESQUITE_MESH_IMPL_DATA_HPP
#define MESQUITE_MESH_IMPL_DATA_HPP

#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "MeshInterface.hpp"

#include <vector>
#include <sys/types.h>

namespace MESQUITE_NS {

class MsqError;

/** Class to store mesh representation for MeshImpl */
class MeshImplData {
  
  public:
  
    MeshImplData() : haveSlavedFlags(false) {}

    /** Clear all data */
    void clear();

    /** Get number of vertices, does not include mid-nodes */
    size_t num_vertices() const;
    /** Get number of elements */
    size_t num_elements() const;
    /** Get number of vertex uses (sum of connectivity length for all elements) 
     * Does not count mid-nodes */
    size_t num_vertex_uses() const;
    
    size_t max_vertex_index() const { return vertexList.size(); }
    size_t max_element_index() const { return elementList.size(); }
    
    /** Copy internal representation into CSR rep 
     *  Does not include mid-nodes. */
    void copy_mesh( size_t* vertex_handle_array,
                    size_t* element_hanlde_array,
                    size_t* element_conn_offsets,
                    size_t* element_conn_indices );
    
    /** Get all vertices, including mid-nodes */                
    void all_vertices( std::vector<size_t>& list, MsqError& err ) const;
    
    /** Get all elements */
    void all_elements( std::vector<size_t>& list, MsqError& err ) const;
    
    /** Check if passed vertex index is valid */
    inline bool is_vertex_valid( size_t index ) const
    {return index < vertexList.size() && vertexList[index].valid; }
    
    /** Check if passed element index is valid */
    inline bool is_element_valid( size_t index ) const
      { return index < elementList.size() &&
               !elementList[index].connectivity.empty(); }

    /** Check if the specified node is used as a mid-node on any element */
    bool is_mid_node( size_t index ) const;
    /** Check if the specified node is used as a corner vertex on any element */
    bool is_corner_node( size_t index ) const;

    /** Get vertex coordinates */
    const Vector3D& get_vertex_coords( size_t index, MsqError& err ) const;
    
    /** Set vertex coordinates */
    void set_vertex_coords( size_t index, 
                            const Vector3D& coords,
                            MsqError& err );
    
    /** Get vertex fixed flag */
    bool vertex_is_fixed( size_t index, MsqError& err ) const;
    
    /** Get vertex slaved flag */
    bool vertex_is_slaved( size_t index, MsqError& err ) const;
    
    /** Set vertex fixed flag */
    void fix_vertex( size_t index, bool flag, MsqError& err );
    
    /** Set vertex slaved flag */
    void slave_vertex( size_t index, bool flag, MsqError& err );
    
    /** Get vertex byte */
    unsigned char get_vertex_byte( size_t index, MsqError& err ) const;
    
    /** Set vertex byte */
    void set_vertex_byte( size_t index, unsigned char value, MsqError& err );
    
    /** Get element type */
    EntityTopology element_topology( size_t index, MsqError& err ) const;
    
    /** Set element type */
    void element_topology( size_t index, EntityTopology type, MsqError& err );
    
    /** Get element connectivity list, including mid-nodes */
    const std::vector<size_t>& element_connectivity( size_t index, MsqError& err ) const;
    
    /** Get vertex adjacency list */
    const std::vector<size_t>& vertex_adjacencies( size_t index, MsqError& err ) const;
    
    /** Allocate space for specified number of vertices */
    void allocate_vertices( size_t count, MsqError& err );
    
    /** Allocate space for specified number of elements */
    void allocate_elements( size_t count, MsqError& err );
    
    /** Set allocated but unset veretx to specified values */
    void reset_vertex( size_t index, 
                       const Vector3D& coords, 
                       bool fixed, 
                       MsqError& err );
    
    /** 
     *  Clear element at specified index (if any) including 
     *  connectivity and adjacency data, and re-initialize with
     *  passed data. 
     */
    void reset_element( size_t index, 
                        const std::vector<long>& vertices,
                        EntityTopology topology,
                        MsqError& err  );
    void reset_element( size_t index, 
                        const std::vector<size_t>& vertices,
                        EntityTopology topology,
                        MsqError& err  );
    
      /** Add a new vertex */
    size_t add_vertex( const Vector3D& coords, bool fixed, MsqError& err );
      /** Add a new element */
    size_t add_element( const std::vector<long>& vertices,
                        EntityTopology topology,
                        MsqError& err  );
    size_t add_element( const std::vector<size_t>& vertices,
                        EntityTopology topology,
                        MsqError& err  );
    
      /** Delete a vertex - may not be referenced by any element */
    void delete_vertex( size_t index, MsqError& err );
      /** Delete an element */
    void delete_element( size_t index, MsqError& err );
    
    /** Get all mid-nodes and their adjacent corner vertices */
    void copy_higher_order( std::vector<size_t>& mid_nodes,
                            std::vector<size_t>& vertices,
                            std::vector<size_t>& vertex_indices,
                            std::vector<size_t>& index_offsets,
                            MsqError& err );
    
    /** \brief Get elements adjacent to ALL of the passed nodes.
     *
     * Return the list of elements that is the intersection of the
     * adjacency lists of the specified vertices.
     */
    void get_adjacent_elements( std::vector<size_t>::const_iterator nodes,
                                std::vector<size_t>::const_iterator nodes_end,
                                std::vector<size_t>& elems_out,
                                MsqError& err );
    
    /**\brief Skin mesh
     *
     * Get the boundary of a mesh as element sides
     *
     *\param sides Element sides as pairs of values : { elem_index, side_number }
     */
    void skin( std::vector<size_t>& sides, MsqError& err );
    
    bool have_slaved_flags() const
      { return haveSlavedFlags; }
  private:
  
    /**\brief helper function for skinning
     *
     * Check if any elements adjacent to a side of an element
     * are of the same dimension as the input element.
     *\param elem  The element 
     *\param nodes The nodes composing the side of the element
     */
    bool has_adjacent_elements( size_t elem,
                                const std::vector<size_t>& nodes,
                                MsqError& err );
     
  
      /** Clear existing element data */
    void clear_element( size_t index, MsqError& err );
    
      /** Set cleared element */
    void set_element( size_t index, 
                      const std::vector<long>& vertices,
                      EntityTopology topology,
                      MsqError& err  );
    
      /** Set cleared element */
    void set_element( size_t index, 
                      const std::vector<size_t>& vertices,
                      EntityTopology topology,
                      MsqError& err  );
  
      /** Struct holding a vertex */
    struct Vertex {
      Vertex( const Vector3D& pos, bool is_fixed )
        : coords(pos), midcount(0), fixed(is_fixed), valid(true), byte('\0') {}
        
      Vertex() : midcount(0), valid(false), byte('\0') {}
    
      Vector3D coords;                     /**< location */
      std::vector<size_t> adjacencies; /**< indices of adjacent elements */
      unsigned midcount;                   /**< num elements referencing this as a mid-node */
      bool fixed;                          /**< is fixed */
      bool slaved;
      bool valid;                          /**< is a valid (initialized) array entry */
      unsigned char byte;                  /**< mark */
    };
    
      /** Struct holding an element */
    struct Element {
      std::vector<size_t> connectivity; /**< list of vertex indices */
      EntityTopology topology;             /**< element type */
      Element()
          : topology(MIXED)
        {}
      
    };
    
    std::vector<Vertex> vertexList;     /**< Array of vertices */
    std::vector<Element> elementList;   /**< Array of elements */
    
    /** List of unused indices in vertex list */
    std::vector<size_t> deletedVertexList;
    /** List of unused indices in element list */
    std::vector<size_t> deletedElementList;
    
    bool haveSlavedFlags;
};

/**\brief VertexIterator for MeshImpl
 *
 * Iterate over valid vertex indices
 */
class MeshImplVertIter : public VertexIterator
{
  private:
    MeshImplData* mesh;
    size_t index;
    
  public:
    
    MeshImplVertIter( MeshImplData* data )
      : mesh(data) { restart(); }
    
    virtual ~MeshImplVertIter();
    
    virtual void restart();
    
    virtual void operator++();
    
    virtual Mesh::VertexHandle operator*() const;
    
    virtual bool is_at_end() const;
    
};


/**\brief ElementIterator for MeshImpl
 *
 * Iterate over valid element indices
 */
class MeshImplElemIter : public ElementIterator
{
  private:
    MeshImplData* mesh;
    size_t index;
    
  public:
    
    MeshImplElemIter( MeshImplData* data )
      : mesh(data) { restart(); }
    
    virtual ~MeshImplElemIter();
    
    virtual void restart();
    
    virtual void operator++();
    
    virtual Mesh::ElementHandle operator*() const;
    
    virtual bool is_at_end() const;
    
};



} // namespace Mesquite

#endif 
