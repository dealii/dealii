// test comment, remove when done
/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file MappingFunction.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "MappingFunction.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "IdealElements.hpp"

namespace MESQUITE_NS {

NodeSet
MappingFunction::sample_points( NodeSet higher_order ) const
{
  higher_order.set_all_corner_nodes( element_topology() );
  return higher_order;
}


void
MappingFunction::convert_connectivity_indices_impl( EntityTopology topo,
                                               int input_type,
                                               int output_type,
                                               size_t* index_list,
                                               unsigned num_indices,
                                               MsqError& err )
{
  bool in_edges, in_faces, in_region, out_edges, out_faces, out_region;
  TopologyInfo::higher_order( topo, input_type, in_edges, in_faces, 
                              in_region, err ); MSQ_ERRRTN(err);
  TopologyInfo::higher_order( topo, output_type, out_edges, out_faces, 
                              out_region, err ); MSQ_ERRRTN(err);

    // We could probably use TopologyInfo to do this more forward-compatible,
    // but for efficiency assume the current ITAPS node ordering, where
    // all mid-edge nodes occur before mid-face nodes and the mid-region
    // node is always last.
    
    // If both have mid-region nodes and they don't have the same stuff 
    // preceeding the mid-region node, then we need to change that index.
  bool region_diff = in_region && out_region && (in_faces != out_faces || in_edges != out_edges);
    // If both have mid-face nodes and one has mid-edge nodes and the other
    // does not, then we need to change the face indices.
  bool face_diff = in_faces && out_faces && in_edges != out_edges;
    // if nothing to change, return
  if (!face_diff && !region_diff)
    return;

  const unsigned corners = TopologyInfo::corners(topo);
  const unsigned edges   = TopologyInfo::edges(topo);
  const unsigned faces   = TopologyInfo::faces(topo);
  const unsigned in_face_offset = in_edges ? corners+edges : corners;
  const unsigned in_regn_offset = in_faces ? in_face_offset+faces : in_face_offset;
  const unsigned out_face_offset = out_edges ? corners+edges : corners;
  const unsigned out_regn_offset = out_faces ? out_face_offset+faces : out_face_offset;
  
    // In the code below, assertions are used to validate the input
    // connectivity data as we assume it is an internal mesquite coding
    // error for it to be inconsistent.  True error checking is used
    // if the elements are incompatible (the index list for the input
    // type contains indices for which there is no correpsonding logical
    // node location in the connectivity list of the output element type)
    // because that indicates an invalid setup (the combination of element
    // type and slave nodes does not result in a reduced element that is
    // compatible with the mapping function.)  The latter should probably
    // have been caught by the mapping function, but to be safe we check
    // again here.
  
  for (size_t i = 0; i < num_indices; ++i) {
    if (index_list[i] < in_face_offset) { // corner or mid-edge node
      // nothing to change for these, but check that if it is a mid-edge
      // node that the other type also has edges
      if (index_list[i] >= corners && !out_edges) {
        MSQ_SETERR(err)("Incompatible nodes present.", MsqError::UNSUPPORTED_ELEMENT );
        return;
      }
    }
    else if (index_list[i] < in_regn_offset) { // mid-face node
      assert( TopologyInfo::dimension(topo) == 3 || index_list[i] == (unsigned)input_type - 1 );
      if (!out_faces) {
        MSQ_SETERR(err)("Incompatible nodes present.", MsqError::UNSUPPORTED_ELEMENT );
        return;
      }
        // working with unsigned type (size_t), so make sure we express this
        // such that there are no intermediate negative values.
      index_list[i] = index_list[i] + out_face_offset - in_face_offset;
    }
    else { // region
      assert( in_region );
      assert( TopologyInfo::dimension(topo) == 3 && index_list[i] == (unsigned)input_type - 1 );
      if (!out_region) {
        MSQ_SETERR(err)("Incompatible nodes present.", MsqError::UNSUPPORTED_ELEMENT );
        return;
      }
        // working with unsigned type (size_t), so make sure we express this
        // such that there are no intermediate negative values.
     index_list[i] = index_list[i] + out_regn_offset - in_regn_offset;
    }
  }
}

void MappingFunction2D::jacobian( const PatchData& pd,
                                  size_t element_number,
                                  NodeSet nodeset,
                                  Sample location,
                                  size_t* vertex_patch_indices_out,
                                  MsqVector<2>* d_coeff_d_xi_out,
                                  size_t& num_vtx_out,
                                  MsqMatrix<3,2>& jacobian_out,
                                  MsqError& err ) const
{
  const MsqMeshEntity& elem = pd.element_by_index( element_number );
  const size_t* conn = elem.get_vertex_index_array();
  
  derivatives( location, nodeset, vertex_patch_indices_out,
               d_coeff_d_xi_out, num_vtx_out, err ); MSQ_ERRRTN(err);
 
  convert_connectivity_indices( elem.node_count(), vertex_patch_indices_out, 
                                num_vtx_out, err );  MSQ_ERRRTN(err);
 
  jacobian_out.zero();
  size_t w = 0;
  for (size_t r = 0; r < num_vtx_out; ++r) {
    size_t i = conn[vertex_patch_indices_out[r]];
    MsqMatrix<3,1> coords( pd.vertex_by_index( i ).to_array() );
    jacobian_out += coords * transpose(d_coeff_d_xi_out[r]);
    if (i < pd.num_free_vertices()) {
      vertex_patch_indices_out[w] = i;
      d_coeff_d_xi_out[w] = d_coeff_d_xi_out[r];
      ++w;
    }
  }
  num_vtx_out = w;
}

void MappingFunction2D::ideal( Sample location,
                               MsqMatrix<3,2>& J,
                               MsqError& err ) const
{
  const Vector3D* coords = unit_element( element_topology(), true );
  if (!coords) {
    MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
    return;
  }
  
  const unsigned MAX_VERTS = 4;
  MsqVector<2> d_coeff_d_xi[MAX_VERTS];
  size_t indices[MAX_VERTS], num_vtx = 0;
  derivatives( location, NodeSet(), indices,
               d_coeff_d_xi, num_vtx, err ); MSQ_ERRRTN(err);
  assert(num_vtx > 0 && num_vtx <= MAX_VERTS);
  
  J.zero();
  for (size_t r = 0; r < num_vtx; ++r) {
    MsqMatrix<3,1> c( coords[indices[r]].to_array() );
    J += c * transpose(d_coeff_d_xi[r]);
  }
  
  double size = sqrt(sqrt(fabs(det(transpose(J) * J))));
  assert(size > -1e-15); // no negative jacobians for ideal elements!
  divide( 1.0, size, size );
  J *= size;
}

void MappingFunction3D::jacobian( const PatchData& pd,
                                  size_t element_number,
                                  NodeSet nodeset,
                                  Sample location,
                                  size_t* vertex_patch_indices_out,
                                  MsqVector<3>* d_coeff_d_xi_out,
                                  size_t& num_vtx_out,
                                  MsqMatrix<3,3>& jacobian_out,
                                  MsqError& err ) const
{
  const MsqMeshEntity& elem = pd.element_by_index( element_number );
  const size_t* conn = elem.get_vertex_index_array();
  
  derivatives( location, nodeset, vertex_patch_indices_out,
               d_coeff_d_xi_out, num_vtx_out, err ); MSQ_ERRRTN(err);
 
  convert_connectivity_indices( elem.node_count(), vertex_patch_indices_out, 
                                num_vtx_out, err );  MSQ_ERRRTN(err);
 
  jacobian_out.zero();
  size_t w = 0;
  for (size_t r = 0; r < num_vtx_out; ++r) {
    size_t i = conn[vertex_patch_indices_out[r]];
    MsqMatrix<3,1> coords( pd.vertex_by_index( i ).to_array() );
    jacobian_out += coords * transpose(d_coeff_d_xi_out[r]);
    if (i < pd.num_free_vertices()) {
      vertex_patch_indices_out[w] = i;
      d_coeff_d_xi_out[w] = d_coeff_d_xi_out[r];
      ++w;
    }
  }
  num_vtx_out = w;
}


void MappingFunction3D::ideal( Sample location,
                               MsqMatrix<3,3>& J,
                               MsqError& err ) const
{
  const Vector3D* coords = unit_element( element_topology(), true );
  if (!coords) {
     MSQ_SETERR(err)(MsqError::UNSUPPORTED_ELEMENT);
     return;
  }

  const unsigned MAX_VERTS = 8;
  MsqVector<3> d_coeff_d_xi[MAX_VERTS];
  size_t indices[MAX_VERTS], num_vtx = 0;
  derivatives( location, NodeSet(), indices,
               d_coeff_d_xi, num_vtx, err ); MSQ_ERRRTN(err);
  assert(num_vtx > 0 && num_vtx <= MAX_VERTS);
  
  J.zero();
  for (size_t r = 0; r < num_vtx; ++r) {
    MsqMatrix<3,1> c( coords[indices[r]].to_array() );
    J += c * transpose(d_coeff_d_xi[r]);
  }
  
  double size = Mesquite::cbrt(fabs(det(J)));
  assert(size > -1e-15); // no negative jacobians for ideal elements!
  divide( 1.0, size, size );
  J *= size;
}

} // namespace Mesquite
