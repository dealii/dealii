/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file EdgeQM.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_EDGE_QM_HPP
#define MSQ_EDGE_QM_HPP

#include "Mesquite.hpp"
#include "QualityMetric.hpp"

namespace MESQUITE_NS {

/**\brief Base type for quality metrics evaluated for each edge */
class EdgeQM : public QualityMetric
{
public:
  
  MESQUITE_EXPORT virtual ~EdgeQM();
  
  MESQUITE_EXPORT virtual MetricType get_metric_type() const
    { return VERTEX_BASED; }
  
  /**\brief Returns list of edge indices in PatchData
   *
   *This method returns metric evaluation points for every
   *logical edge in the patch if \c free_vertices_only is
   *false.  If \c free_vertices_only is true then only the 
   *subset of edges adjacent to at least one free vertex are
   *returned.
   */
  MESQUITE_EXPORT virtual 
  void get_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles,
                        bool free_vertices_only, 
                        MsqError& err );
  
  /**\brief Returns list of edge indices in PatchData 
   *
   *This method returns metric evaluation points only a subset
   *of the logical edges in a patch such that if one iterates 
   *over the mesh using element-on-vertex patches a given edge
   *is returned only once for the set of all patches.  This is
   *accomplished by returning only edges adjacent to vertices
   *without the MSQ_PATCH_FIXED flag set, and only if the handle for
   *the opposite vertex is greater than the one with the flag
   *set.
   */
  MESQUITE_EXPORT virtual 
  void get_single_pass( PatchData& pd, 
                        std::vector<size_t>& handles,
                        bool free_vertices_only, 
                        MsqError& err );

  MESQUITE_EXPORT static
  void get_edge_evaluations( PatchData& pd, 
                             std::vector<size_t>& handles,
                             bool free_vertices_only, 
                             bool single_pass_evaluate,
                             MsqError& err );

   /**\brief Default implementation for all edge-based metrics
    *
    * Fill 'indices' with all free vertex indices in element,
    * and call 'evaluate'.
    */
  MESQUITE_EXPORT virtual
  bool evaluate_with_indices( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 MsqError& err );
                 
  enum {
    ELEM_EDGE_BITS = 4,
    ELEM_INDEX_BITS = 8*sizeof(size_t) - ELEM_EDGE_BITS,
    ELEM_EDGE_MASK = (((size_t)1) << ELEM_INDEX_BITS) - 1
  };
  
  inline static size_t handle( unsigned edge_no, size_t elem_idx )
    { assert(elem_idx <= ELEM_EDGE_MASK);
      return (((size_t)edge_no) << ELEM_INDEX_BITS) | elem_idx; }
      
  inline static unsigned edge( size_t handle )
    { return handle >> ELEM_INDEX_BITS; }
  
  inline static unsigned elem( size_t handle )
    { return handle & ELEM_EDGE_MASK; }
};



} // namespace MESQUITE_NS

#endif
