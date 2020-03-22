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

#ifndef MESQUITE_TOPOLOGY_INFO_HPP
#define MESQUITE_TOPOLOGY_INFO_HPP

#include "Mesquite.hpp"
#include "Sample.hpp"
#include <string.h>

namespace MESQUITE_NS
{

class MsqError;

/** \brief Information about different element topologies */
class MESQUITE_EXPORT TopologyInfo
{
  public:
  
    static const char* name( EntityTopology topo )
      { return topo > MIXED ? 0 : instance.longNames[topo]; }
    static const char* short_name( EntityTopology topo )
      { return topo > MIXED ? 0 : instance.shortNames[topo]; }
    
      /** \brief Dimension of element topology */
    static unsigned dimension( EntityTopology topo )
      { return topo >= MIXED ? 0: instance.dimMap[topo]; }
    
      /** \brief Number of adjacent entities of a specified dimension 
       *
       * For a given element topology, get the number of adjacent entities
       * of the specified dimension.  
       */
    static unsigned adjacent( EntityTopology topo, unsigned dimension )
      { return (topo >= MIXED) ? 0 : instance.adjMap[topo][dimension]; }
    
      /** \brief Get number of sides a given topology type has 
       *
       * Get the number of sides for a given element topology.  Returns
       * the number of adjacent entities of one less dimension.  The number
       * of faces for a volume element and the number of edges for a face
       * element 
       */
    static unsigned sides( EntityTopology topo )
      { return (topo >= MIXED || instance.dimMap[topo] < 1) ? 0 : 
          instance.adjMap[topo][instance.dimMap[topo]-1]; }
    
      /** \brief Get the number of defining vertices for a given element topology 
       *
       * Get the number of corner vertices necessary to define an element
       * of the specified topology.  This is the number of nodes a linear
       * element of the specified topology will have.
       */
    static unsigned corners( EntityTopology topo )
      { return adjacent(topo, 0); }
    
      /** \brief Get the number of edges in a given topology */
    static unsigned edges( EntityTopology topo )
      { return adjacent(topo, 1); }
      
      /** \brief Get the number of faces in a given topology */
    static unsigned faces( EntityTopology topo )
      { return adjacent(topo, 2 ); }
      
      /** \brief Check which mid-nodes a higher-order element has.
       *
       * Assuming at most one mid-node per sub-entity per dimension
       * (i.e. at most one mid-node per edge, at most one mid-node per face, etc.)
       * determine which mid-nodes are present given the topology
       * and total number of nodes.
       */
    static void higher_order( EntityTopology topo, unsigned num_nodes,
                              bool& midedge, bool& midface, bool& midvol,
                              MsqError &err );
      
      /** \brief Check which mid-nodes a higher-order element has.
       *
       * Assuming at most one mid-node per sub-entity per dimension
       * (i.e. at most one mid-node per edge, at most one mid-node per face, etc.)
       * determine which mid-nodes are present given the topology
       * and total number of nodes.  This function is similar to the
       * previous one, except that that it returns a set of bits, one per
       * side dimension, rather than separate bool values.  If the bit at position
       * one (the second least significant bit) has a value of one, then the
       * element has mid-edge nodes.  If the bit at position two (the third to
       * least signficiant bit) has a value of one then the element has mid-face
       * nodes.
       *\code
       *  int ho = TopologyInfo::higher_order( topo, num_nodes, err );
       *  bool have_midedge = !!(ho & 1<<1);
       *  bool have_midface = !!(ho & 1<<2);
       *  bool have_midvol  = !!(ho & 1<<3);
       *\endocde
       *
       * The advantange of this form of the function over the previous is 
       * that a) it is possible to check for mid-nodes on sub-entities of
       * a varialbe dimension 'd':
       *\code
       *  if (ho & (1<<d)) { ... }
       *\code
       * and b) it is convienent to test if an element has any higher-order
       * nodes:
       *\code
       *  int ho = TopologyInfo::higher_order( topo, num_nodes, err );
       *  if (!ho) // if linear element
       *    { ... }
       *\endocde        
       */
    static int higher_order( EntityTopology topo, unsigned num_nodes, MsqError &err );
    
      /**\brief Given a side, return index of mid-vertex for that side.
       *
       * Given a side specification (e.g. the first edge), return the
       * index of of the correponding mid-side node in the canoncial
       * ordering of the element connectivity.  Returns -1 if the element
       * doesn't have the specified mid-side node. 
       *
       *\param topo   The element topology
       *\param num_nodes  The number of nodes in the element type.
       *\param side_dimension  The dimension of the side (e.g. 1 = edge, 2 = face)
       *\param side_number     The number of the side (e.g. 0 for first edge/face, etc.)
       *\return  Index (zero-based position) of higher-order node in canonical 
       *         ordering of element connectivity, or -1 of element type contains
       *         no such node.
       */
    static int higher_order_from_side( EntityTopology topo,
                                       unsigned num_nodes,
                                       unsigned side_dimension,
                                       unsigned side_number,
                                       MsqError& err );
    
      /**\brief Get side given a higher-order node */
    static void side_from_higher_order( EntityTopology topo,
                                        unsigned num_nodes,
                                        unsigned node_number,
                                        unsigned& side_dim_out,
                                        unsigned& side_num_out,
                                        MsqError& err );

    /** Get logical position given an element type node node index*/
    static inline
    Sample sample_from_node( EntityTopology topo,
                             unsigned num_nodes,
                             unsigned node_number,
                             MsqError& err )
    {
      unsigned dim, num;
      side_from_higher_order( topo, num_nodes, node_number, dim, num, err );
      return Sample(dim, num);
    }
    /** Get node index from logical position */
    static inline
    int node_from_sample( EntityTopology topo, 
                          unsigned num_nodes,
                          Sample sample,
                          MsqError& err )
    {
      return higher_order_from_side( topo, num_nodes, sample.dimension,
                                     sample.number, err );
    }
      /**\brief Get indices of edge ends in element connectivity array 
       *
       * Given an edge number in (0,edges(type)], return which positions
       * in the connectivity list for the element type correspond to the
       * end vertices of that edge.
       */
    static const unsigned* edge_vertices( EntityTopology topo,
                                          unsigned edge_number,
                                          MsqError& err );
    static const unsigned* edge_vertices( EntityTopology topo,
                                          unsigned edge_number );
    
     /**\brief Get face corner indices in element connectivity array 
       *
       * Given an face number in (0,faces(type)], return which positions
       * in the connectivity list for the element type correspond to the
       * vertices of that face, ordered in a counter-clockwise cycle
       * around a vector pointing out of the element for an ideal element.
       */
    static const unsigned* face_vertices( EntityTopology topo,
                                          unsigned face_number,
                                          unsigned& num_vertices_out,
                                          MsqError& err );
    static const unsigned* face_vertices( EntityTopology topo,
                                          unsigned face_number,
                                          unsigned& num_vertices_out );

    /**\brief Get corner indices of side 
     *
     * Get the indices into element connectivity list for the 
     * corners/ends of the specified side of the element.  
     * edge_vertices() and face_vertices() are special cases
     * of this method.  
     *
     * If the passed dimension equals that of the specified topology,
     * the side number is ignored and all the corners of the 
     * element are returned.  Fails if side dimension
     * greater than the dimension of the specified topology type.
     */
    static const unsigned* side_vertices( EntityTopology topo,
                                          unsigned side_dimension, 
                                          unsigned side_number,
                                          unsigned& num_verts_out,
                                          MsqError& err );
    static const unsigned* side_vertices( EntityTopology topo,
                                          unsigned side_dimension, 
                                          unsigned side_number,
                                          unsigned& num_verts_out );

     
      
      /**\brief Return which side the specified mid-node lies on 
       *
       * Given an non-linear element type (specified by the
       * topology and length of the connectiivty array) and the 
       * index of a node in the element's connectivity array,
       * return the lower-dimension entity (side) of the element
       * the mid-node lies on.
       *
       *\param topo  Element topology
       *\param connectivity_length Number of nodes in element
       *\param node_index Which node of the element
       *\param side_dimension_out The dimension of the side containing the
       *             midnode (0 = vertex, 1 = edge, 2 = face, 3 = volume)
       *\param side_number_out The canonical number of the side 
       */
    static void side_number( EntityTopology topo, 
                             unsigned connectivity_length,
                             unsigned node_index,
                             unsigned& side_dimension_out,
                             unsigned& side_number_out,
                             MsqError& err );

    /**\brief  Get adjacent corner vertices
     *
     * Given the index of a vertex in an element, get the list of 
     * indices corresponding to the adjacent corner vertices.
     *
     * Adjcent corner vertex indices are returned in the proper
     * order for constructing the active matrix for the corner.
     *
     * Given the array v of all vertices in the patch, the array v_i
     * containing the connectivity list for an element as
     * indices into v, and adj as the result of this function for some
     * corner of the element, the corresponding active matrix A for
     * that corner can be constructed as:
     *  Matrix3D A;
     *  A.set_column( 0, v[v_i[adj[0]]] - v[v_i[0]] );
     *  A.set_column( 1, v[v_i[adj[1]]] - v[v_i[0]] );
     *  A.set_column( 2, v[v_i[adj[2]]] - v[v_i[0]] );
     *
     *\param topo  The element type
     *\param index The index of a corner vertex
     *\param num_adj_out The number of adjacent vertices (output)
     *\return The array of vertex indices
     */
    static const unsigned* adjacent_vertices( EntityTopology topo,
                                              unsigned index,
                                              unsigned& num_adj_out );
     
    /**\brief  Get reverse adjacency offsets
     *
     * Get reverse mapping of results from adjacent_vertices().
     *
     * Let i be the input vertex index.  For each vertex index j
     * for which the result of adjacent_vertices() contains i, return
     * the offset into that result at which i would occur.  The
     * results are returned in the same order as each j is returned
     * in the results of adjacent_vertices(...,i,...).  Thus the
     * combination of the results of adjacent_vertices(...,i,...)
     * and this method provide a reverse mapping of the results of
     * adjacent_vertices(...,j,...) for i in all j.
     * 
     * Given:
     *   const unsigned *a, *b, *r;
     *   unsigned n, nn, c = corners(type);
     *   a = adjacent_vertices( type, i, n );            // for any i < c
     *   r = reverse_vertex_adjacency_offsets( type, i, n );
     *   b = adjacent_vertices( type, a[k], nn );        // for any k < n
     * Then:
     *   b[r[k]] == i
     */
    static const unsigned* reverse_vertex_adjacency_offsets( 
                                              EntityTopology topo,
                                              unsigned index,
                                              unsigned& num_idx_out );
    
      /**\brief Find which edge of an element has the passed vertex indices
      *
      * Find which edge of the element cooresponds to a list of positions
      * in the canonical element ordering.
      *\param topo            The element type
      *\param edge_vertices   The array of side vertex indices
      *\param reversed_out    True if edge is reversed wrt edge_vertices
      *\return                The edge number.
      */
    static unsigned find_edge( EntityTopology topo,
                               const unsigned* edge_vertices,
                               bool& reversed_out,
                               MsqError& err );
    
      /**\brief Find which face of an element has the passed vertex indices
       *
       * Find which face of the element cooresponds to a list of positions
       * in the canonical element ordering.
       *\param topo           The element type
       *\param face_vertices  The array of face vertex indices
       *\param num_face_vertices   The length of face_vertices
       *\param reversed_out   True if face is reversed wrt face_vertices
       *\return               The face number.
       */
    static unsigned find_face( EntityTopology topo,
                               const unsigned* face_vertices,
                               unsigned num_face_vertices,
                               bool& reversed_out,
                               MsqError& err );
  
      /**\brief Find which side of an element has the passed vertex indices
       *
       * Find which side of the element cooresponds to a list of positions
       * in the canonical element ordering.
       *\param topo           The element type
       *\param side_vertices  The array of side vertex indices
       *\param num_vertices   The length of side_vertices
       *\param dimension_out  The dimension of the side
       *\param number_out     The enumerated index for the side
       *\param reversed_out   True if side is reversed wrt side_vertices
       */
    static void find_side( EntityTopology topo, 
                           const unsigned* side_vertices,
                           unsigned num_vertices,
                           unsigned& dimension_out,
                           unsigned& number_out,
                           bool& reversed_out,
                           MsqError& err );

      /**\brief Test if two elements share lower-order topology
       *
       * Test if two elements share lower-order topology (e.g.
       * whether or not two tetrahedra share an edge.)
       *
       * That is compare the 'element_1_side_number'-th lower order 
       * topology of dimension 'side_dimension' on element 1 with the 
       * 'element_2_side_number'-th lower order topology of dimension 
       *'side_dimension' on element 2
       *
       *\param element_1_vertices    The connectivity of the first element
       *\param element_1_topology    The type of the first element
       *\param element_1_side_number Which lower-order topology to compare
       *\param element_2_vertices    The connectivity of the second element
       *\param element_2_topology    The type of the second element
       *\param element_2_side_number Whcih lower-order topology to compare
       *\param side_dimension        The dimension of the lower-order topology
       */
    static bool compare_sides( const size_t* element_1_vertices,
                               EntityTopology element_1_topology,
                               unsigned element_1_side_number,
                               const size_t* element_2_vertices,
                               EntityTopology element_2_topology,
                               unsigned element_2_side_number,
                               unsigned side_dimension,
                               MsqError& err );

  private:
 
    enum {
      MAX_CORNER = 8,
      MAX_EDGES = 12,
      MAX_FACES = 6,
      MAX_FACE_CONN = 5,
      MAX_VERT_ADJ = 4,
      FIRST_FACE = TRIANGLE,
      LAST_FACE = QUADRILATERAL,
      FIRST_VOL= TETRAHEDRON,
      LAST_VOL = PYRAMID
    };
    
    unsigned char dimMap[MIXED];    /**< Get dimension of entity given topology */
    unsigned char adjMap[MIXED][4]; /**< Get number of adj entities of dimension 0, 1 and dimension 2 */
    /** Vertex indices for element edges */
    unsigned edgeMap[LAST_VOL-FIRST_FACE+1][MAX_EDGES][2] ;
    /** Vertex indices for element faces */
    unsigned faceMap[LAST_VOL-FIRST_VOL+1][MAX_FACES][MAX_FACE_CONN];
    /** Vertex-Vertex adjacency map */
    unsigned vertAdjMap[LAST_VOL-FIRST_FACE+1][MAX_CORNER][MAX_VERT_ADJ+1];
    /** Reverse Vertex-Vertex adjacency index map */
    unsigned revVertAdjIdx[LAST_VOL-FIRST_FACE+1][MAX_CORNER][MAX_VERT_ADJ+1];

    const char* longNames[MIXED+1];
    const char* shortNames[MIXED+1];

    TopologyInfo();
    
    static TopologyInfo instance;
  
};

} //namespace Mesquite

#endif
