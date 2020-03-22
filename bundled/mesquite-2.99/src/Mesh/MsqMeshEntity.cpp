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
//
// ORIG-DATE: 16-May-02 at 10:26:21
//  LAST-MOD: 18-Jun-04 at 11:36:07 by Thomas Leurent
//
/*! \file MsqMeshEntity.cpp

\brief All elements in Mesquite are of type MsqMeshEntity. Their associated
functionality is implemented in this file. 
  
    \author Thomas Leurent
    \author Michael Brewer
    \author Darryl Melander
    \date 2002-05-16  
 */

#include "Mesquite.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"

#include <vector>
#include <ostream>
using std::vector;
using std::ostream;
using std::endl;

namespace MESQUITE_NS {

//! Gets the indices of the vertices of this element.
//! The indices are only valid in the PatchData from which
//! this element was retrieved.
//! The order of the vertices is the canonical order for this
//! element's type.
void MsqMeshEntity::get_vertex_indices(std::vector<std::size_t> &vertices) const
{
  vertices.resize( vertex_count() );
  std::copy( vertexIndices, vertexIndices + vertex_count(), vertices.begin() );
}

//! Gets the indices of the vertices of this element.
//! The indices are only valid in the PatchData from which
//! this element was retrieved.
//! The order of the vertices is the canonical order for this
//! element's type.
//! The indices are placed appended to the end of the list.
//! The list is not cleared before appending this entity's vertices.
void MsqMeshEntity::append_vertex_indices(std::vector<std::size_t> &vertex_list) const
{
  vertex_list.insert(vertex_list.end(),
                     vertexIndices,
                     vertexIndices + vertex_count());
}

void MsqMeshEntity::get_node_indices( std::vector<std::size_t>& indices ) const
{
  indices.resize( node_count() );
  std::copy( vertexIndices, vertexIndices + node_count(), indices.begin() );
}

void MsqMeshEntity::append_node_indices( std::vector<std::size_t>& indices ) const
{
  indices.insert( indices.end(), vertexIndices, vertexIndices + node_count() );
}


/*! The centroid of an element containing n vertices with equal masses is located at
  \f[ \b{x} = \frac{ \sum_{i=1}^{n} \b{x}_i }{ n }  \f]
  where \f$ \b{x}_i  ,\, i=1,...,n\f$ are the vertices coordinates.
*/
void MsqMeshEntity::get_centroid(Vector3D &centroid, const PatchData &pd, MsqError &err) const
{
  centroid = 0.0;
  const MsqVertex* vtces = pd.get_vertex_array(err); MSQ_ERRRTN(err);
  size_t nve = vertex_count();
  for (size_t i=0; i<nve; ++i)
    centroid += vtces[vertexIndices[i]];
  centroid /= nve;
}
  

static inline double corner_volume( const Vector3D& v0,
                                    const Vector3D& v1,
                                    const Vector3D& v2,
                                    const Vector3D& v3 )
{
  return (v1 - v0) * (v2 - v0) % (v3 - v0);
}

/*!
  \brief Computes the area of the given element.  Returned value is
  always non-negative.  If the entity passed is not a two-dimensional
  element, an error is set.*/
double MsqMeshEntity::compute_unsigned_area(PatchData &pd, MsqError &err) {
  const MsqVertex* verts=pd.get_vertex_array(err);MSQ_ERRZERO(err);
  double tem=0.0;
  switch (mType)
  {
   
    case TRIANGLE:
      tem =  ((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
              (verts[vertexIndices[2]]-verts[vertexIndices[0]])).length();
      return 0.5*tem;
      
    case QUADRILATERAL:
      tem = ((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
             (verts[vertexIndices[3]]-verts[vertexIndices[0]])).length();
      tem += ((verts[vertexIndices[3]]-verts[vertexIndices[2]])*
              (verts[vertexIndices[1]]-verts[vertexIndices[2]])).length();
      return (tem/2.0);
    
    case POLYGON: 
        // assume convex
      for (unsigned i = 1; i < numVertexIndices-1; ++i)
        tem += ((verts[vertexIndices[i]] - verts[vertexIndices[0]]) *
                (verts[vertexIndices[i+1]] - verts[vertexIndices[0]])).length();
      return 0.5 * tem;
    
    case TETRAHEDRON:
      return 1.0/6.0 * fabs( corner_volume( verts[vertexIndices[0]],
                                            verts[vertexIndices[1]],
                                            verts[vertexIndices[2]],
                                            verts[vertexIndices[3]] ) );
    
    case PYRAMID: {
      Vector3D m = verts[vertexIndices[0]] + verts[vertexIndices[1]] 
                 + verts[vertexIndices[2]] + verts[vertexIndices[3]];
      Vector3D t1 = verts[vertexIndices[0]] - verts[vertexIndices[2]];
      Vector3D t2 = verts[vertexIndices[1]] - verts[vertexIndices[3]];
      tem = ((t1 + t2) * (t1 - t2)) % (verts[vertexIndices[4]] - 0.25 * m);
      return (1.0/12.0) * fabs(tem);
    }
    
    case PRISM: {
      tem  = corner_volume( verts[vertexIndices[0]],
                            verts[vertexIndices[1]],
                            verts[vertexIndices[2]],
                            verts[vertexIndices[3]] );

      tem += corner_volume( verts[vertexIndices[1]],
                            verts[vertexIndices[2]],
                            verts[vertexIndices[3]],
                            verts[vertexIndices[4]] );

      tem += corner_volume( verts[vertexIndices[2]],
                            verts[vertexIndices[3]],
                            verts[vertexIndices[4]],
                            verts[vertexIndices[5]] );
                        
      return 1.0/6.0 * fabs(tem);
    }
    
    case HEXAHEDRON: {
    
      tem  = corner_volume( verts[vertexIndices[1]],
                            verts[vertexIndices[2]],
                            verts[vertexIndices[0]],
                            verts[vertexIndices[5]] );
                            
      tem += corner_volume( verts[vertexIndices[3]],
                            verts[vertexIndices[0]],
                            verts[vertexIndices[2]],
                            verts[vertexIndices[7]] );
                            
      tem += corner_volume( verts[vertexIndices[4]],
                            verts[vertexIndices[7]],
                            verts[vertexIndices[5]],
                            verts[vertexIndices[0]] );
                            
      tem += corner_volume( verts[vertexIndices[6]],
                            verts[vertexIndices[5]],
                            verts[vertexIndices[7]],
                            verts[vertexIndices[2]] );
                            
      tem += corner_volume( verts[vertexIndices[5]],
                            verts[vertexIndices[2]],
                            verts[vertexIndices[0]],
                            verts[vertexIndices[7]] );

      return (1.0/6.0) * fabs(tem);
    }
    
    default:
      MSQ_SETERR(err)("Invalid type of element passed to compute unsigned area.",
                      MsqError::UNSUPPORTED_ELEMENT);
      return 0;
  }
  return 0;
}
                                            
/*!
  \brief Computes the area of the given element.  Returned value can be
  negative.  If the entity passed is not a two-dimensional element, an
  error is set.*/
double MsqMeshEntity::compute_signed_area(PatchData &pd, MsqError &err) {
  const MsqVertex* verts=pd.get_vertex_array(err);MSQ_ERRZERO(err);
  double tem=0.0;
  double tem2=0.0;
  Vector3D surface_normal;
  Vector3D cross_vec;
  size_t element_index=pd.get_element_index(this);
  
  switch (mType)
  {
    
    case TRIANGLE:
      cross_vec=((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
      (verts[vertexIndices[2]]-verts[vertexIndices[0]]));
      pd.get_domain_normal_at_element(element_index,surface_normal,err);
      MSQ_ERRZERO(err);
      tem =  (cross_vec.length()/2.0);
        //if normals do not point in same general direction, negate area
      if(cross_vec%surface_normal<0){ 
        tem *= -1;
      }
      
      return tem;
      
    case QUADRILATERAL:
      cross_vec=((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
                 (verts[vertexIndices[3]]-verts[vertexIndices[0]]));
      pd.get_domain_normal_at_element(element_index,surface_normal,err);
      MSQ_ERRZERO(err);
      tem =  (cross_vec.length()/2.0);
        //if normals do not point in same general direction, negate area
      if(cross_vec%surface_normal<0){ 
        tem *= -1;
      }
      cross_vec=((verts[vertexIndices[3]]-verts[vertexIndices[2]])*
                 (verts[vertexIndices[1]]-verts[vertexIndices[2]]));
      tem2 =  (cross_vec.length()/2.0);
        //if normals do not point in same general direction, negate area
      if(cross_vec%surface_normal<0){ 
        tem2 *= -1;
          //test to make sure surface normal existed
          //if(surface_normal.length_squared()<.5){
          //err.set_msg("compute_signed_area called without surface_normal available.");
          //}  
      }
      return (tem + tem2);
      
    default:
      MSQ_SETERR(err)("Invalid type of element passed to compute unsigned area.",
                       MsqError::UNSUPPORTED_ELEMENT);
      return 0;
  };
  return 0.0;
}
    
/*!Appends the indices (in the vertex array) of the vertices to connected
  to vertex_array[vertex_index] to the end of the vector vert_indices.
  The connected vertices are right-hand ordered as defined by the
  entity.
  
*/
void MsqMeshEntity::get_connected_vertices( 
                                    std::size_t vertex_index,
                                    std::vector<std::size_t> &vert_indices,
                                    MsqError &err)
{
    //index is set to the index in the vertexIndices corresponding
    //to vertex_index
  int index;
  for (index = vertex_count() - 1; index >= 0; --index)
    if (vertexIndices[index] == vertex_index)
      break;
  if (index < 0)
  {
    MSQ_SETERR(err)("Invalid vertex index.", MsqError::INVALID_ARG);
    return;
  }
  
  unsigned n;
  const unsigned* indices = TopologyInfo::adjacent_vertices( mType, index, n );
  if (!indices)
    MSQ_SETERR(err)("Element type not available", MsqError::INVALID_ARG);
  for (unsigned i = 0; i < n; ++i)
    vert_indices.push_back( vertexIndices[indices[i]] );
}

/*! Gives the normal at the surface point corner_pt ... but if not available,
    gives the normalized cross product of corner_vec1 and corner_vec2. 
  */
/*
void MsqMeshEntity::compute_corner_normal(size_t corner,
                                          Vector3D &normal,
                                          PatchData &pd, 
                                          MsqError &err)
{
  if ( get_element_type()==TRIANGLE || get_element_type()==QUADRILATERAL ) 
  {
      // There are two cases where we cannot get a normal from the 
      // geometry that are not errors:
      // 1) There is no domain set
      // 2) The vertex is at a degenerate point on the geometry (e.g. 
      //     tip of a cone.)
    
    bool have_normal = false;

      // Get normal from domain
    if (pd.domain_set()) 
    {
      size_t index = pd.get_element_index(this);
      pd.get_domain_normal_at_corner( index, corner, normal, err );
      MSQ_ERRRTN(err);
      
      double length = normal.length();
      if (length > DBL_EPSILON) 
      {
        have_normal = true;
        normal /= length;
      }
    }

      // If failed to get normal from domain, calculate
      // from adjacent vertices.
    if ( !have_normal )
    {
      const int num_corner = this->vertex_count();
      const int prev_corner = (corner + num_corner - 1) % num_corner;
      const int next_corner = (corner + 1 ) % num_corner;
      const size_t this_idx = vertexIndices[corner];
      const size_t prev_idx = vertexIndices[prev_corner];
      const size_t next_idx = vertexIndices[next_corner];
      normal = (pd.vertex_by_index(next_idx) - pd.vertex_by_index(this_idx))
             * (pd.vertex_by_index(prev_idx) - pd.vertex_by_index(this_idx));
      normal.normalize();
    }
  }
  else
    MSQ_SETERR(err)("Should only be used for faces (tri, quads, ...).",
                    MsqError::INVALID_ARG);
}
*/

void MsqMeshEntity::compute_corner_normals( Vector3D normals[],
                                            PatchData &pd, 
                                            MsqError &err)
{
  EntityTopology type = get_element_type();
  if (type != TRIANGLE && type != QUADRILATERAL && type != POLYGON)
  {
      MSQ_SETERR(err)("Should only be used for faces (tri, quads, ...).",
                    MsqError::INVALID_ARG);
      return;
  }
  
  
    // There are two cases where we cannot get a normal from the 
    // geometry that are not errors:
    // 1) There is no domain set
    // 2) The vertex is at a degenerate point on the geometry (e.g. 
    //     tip of a cone.)

    // Get normal from domain
  if (pd.domain_set()) 
  {
    size_t index = pd.get_element_index(this);
    pd.get_domain_normals_at_corners( index, normals, err );
    MSQ_ERRRTN(err);
  }

    // Check if normals are valid (none are valid if !pd.domain_set())
  const unsigned count = vertex_count();
  for (unsigned i = 0; i < count; ++i)
  {
      // If got valid normal from domain, 
      // make it a unit vector and continue.
    if (pd.domain_set()) 
    {
      double length = normals[i].length();
      if (length > DBL_EPSILON)
      {
        normals[i] /= length;
        continue;
      }
    }
    
    const size_t prev_idx = vertexIndices[(i + count - 1) % count];
    const size_t this_idx = vertexIndices[i];
    const size_t next_idx = vertexIndices[(i + 1) % count];

      // Calculate normal using edges adjacent to corner
    normals[i] = (pd.vertex_by_index(next_idx) - pd.vertex_by_index(this_idx))
               * (pd.vertex_by_index(prev_idx) - pd.vertex_by_index(this_idx));
    normals[i].normalize();
  }
}

ostream& operator<<( ostream& stream, const MsqMeshEntity& entity )
{
  size_t num_vert = entity.vertex_count();
  stream << "MsqMeshEntity " << &entity << " with vertices ";
  for (size_t i = 0; i < num_vert; ++i)
    stream << entity.vertexIndices[i] << " ";
  stream << endl;
  return stream;
}




/*! Get a array of indices that specifies for the given vertex
  the correct matrix map.  This is used by the I_DFT point
  relaxation methods in the laplacian smoothers.
  
*/
size_t MsqMeshEntity::get_local_matrix_map_about_vertex(
  PatchData &pd, MsqVertex* vert, size_t local_map_size,
  int* local_map, MsqError &err) const
{
    //i iterates through elem's vertices
  int i=0;
    //index is set to the index in the vertexIndices corresponding
    //to vertex_index
  int index=-1;
  int return_val=0;
  const MsqVertex* vertex_array = pd.get_vertex_array(err);
  if(err)
    return return_val;
  
  switch (mType)
  {
    case TRIANGLE:
      MSQ_SETERR(err)("Requested function not yet supported for Triangles.",
                      MsqError::NOT_IMPLEMENTED);
      
      break;
      
    case QUADRILATERAL:
      MSQ_SETERR(err)("Requested function not yet supported for Quadrilaterals.",
                      MsqError::NOT_IMPLEMENTED);
       
      break;
      
    case TETRAHEDRON:
      if(local_map_size<4){
        MSQ_SETERR(err)("Array of incorrect length sent to function.",
                        MsqError::INVALID_ARG);
        return return_val;
      }
      return_val = 4;
      while(i<4)
      {
        if(&vertex_array[vertexIndices[i]]==vert)
        {
          index=i;
          break;
        }
        ++i;
      }
      switch(index){
        case(0):
          local_map[0]=0;
          local_map[1]=1;
          local_map[2]=2;
          local_map[3]=3;
          break;
        case(1):
          local_map[0]=1;
          local_map[1]=0;
          local_map[2]=3;
          local_map[3]=2;
          break;
        case(2):
          local_map[0]=2;
          local_map[1]=3;
          local_map[2]=0;
          local_map[3]=1;
          break;
        case(3):
          local_map[0]=3;
          local_map[1]=2;
          local_map[2]=1;
          local_map[3]=0;
          break;
        default:
          local_map[0]=-1;
          local_map[1]=-1;
          local_map[2]=-1;
          local_map[3]=-1;
      };
      
      break;
      
    case HEXAHEDRON:
      MSQ_SETERR(err)("Requested function not yet supported for Hexahedrons.",
                      MsqError::NOT_IMPLEMENTED);
      
      break;
    default:
      MSQ_SETERR(err)("Element type not available", MsqError::UNSUPPORTED_ELEMENT);
      break;
  }
  return return_val;
  
}


void MsqMeshEntity::check_element_orientation(
  PatchData &pd, int& inverted, int& total, MsqError &err)
{
  NodeSet all = all_nodes( err ); MSQ_ERRRTN(err);
  unsigned i;

  if (TopologyInfo::dimension(mType) == 2) {
    if (!pd.domain_set()) {
      total = 0;
      inverted = 0;
      return;
    }
  
    const MappingFunction2D* mf = pd.get_mapping_function_2D( mType );
    if (!mf) {
      check_element_orientation_corners( pd, inverted, total, err );
      return;
    }
   
    NodeSet sample = mf->sample_points( all );
    total = sample.num_nodes();
    inverted = 0;
    
    if (sample.have_any_corner_node()) {
      for (i = 0; i < TopologyInfo::corners(mType); ++i)
        if (sample.corner_node(i)) 
          inverted += inverted_jacobian_2d( pd, all, Sample(0,i), err );
    }
    if (sample.have_any_mid_edge_node()) {
      for (i = 0; i < TopologyInfo::edges(mType); ++i)
        if (sample.mid_edge_node(i))
          inverted += inverted_jacobian_2d( pd, all, Sample(1,i), err );
    } 
    if (sample.have_any_mid_face_node())
      inverted += inverted_jacobian_2d( pd, all, Sample(2,0), err );
  }
  else {
    const MappingFunction3D* mf = pd.get_mapping_function_3D( mType );
    if (!mf) {
      check_element_orientation_corners( pd, inverted, total, err );
      return;
    }
    
    NodeSet sample = mf->sample_points( all );
    total = sample.num_nodes();
    inverted = 0;
    
    if (sample.have_any_corner_node()) {
      for (i = 0; i < TopologyInfo::corners(mType); ++i)
        if (sample.corner_node(i))
          inverted += inverted_jacobian_3d( pd, all, Sample(0,i), err );
    }
    if (sample.have_any_mid_edge_node()) {
      for (i = 0; i < TopologyInfo::edges(mType); ++i)
        if (sample.mid_edge_node(i))
          inverted += inverted_jacobian_3d( pd, all, Sample(1,i), err );
    } 
    if (sample.have_any_mid_face_node()) {
      for (i = 0; i < TopologyInfo::faces(mType); ++i)
        if (sample.mid_face_node(i))
          inverted += inverted_jacobian_3d( pd, all, Sample(2,i), err );
    } 
    if (sample.have_any_mid_region_node()) {
      inverted += inverted_jacobian_3d( pd, all, Sample(3,0), err );
    }
  }
}

bool
MsqMeshEntity::inverted_jacobian_3d( PatchData& pd, NodeSet nodes, Sample sample, MsqError& err )
{
  MsqMatrix<3,3> J;
  MsqVector<3> junk[27];
  size_t junk2[27], junk3;
  assert(node_count() <= 27);

  const MappingFunction3D* mf = pd.get_mapping_function_3D( mType );
  mf->jacobian( pd, pd.get_element_index(this), nodes, 
                sample, junk2, junk, junk3, J, err );
  MSQ_ERRZERO(err);
  //const double size_eps_sqr = sqr_Frobenius( J ) * DBL_EPSILON;
  const double d = det(J);
  double l1 = J.column(0) % J.column(0);
  double l2 = J.column(1) % J.column(1);
  double l3 = J.column(2) % J.column(2);
  return d < 0 || d*d < DBL_EPSILON*DBL_EPSILON * l1*l2*l3;
}

bool 
MsqMeshEntity::inverted_jacobian_2d( PatchData& pd, NodeSet nodes, Sample sample, MsqError& err )
{
  MsqMatrix<3,2> J;
  MsqVector<2> junk[9];
  size_t junk2[9], junk3;
  assert(node_count() <= 9);

  const int idx = pd.get_element_index(this);
  const MappingFunction2D* mf = pd.get_mapping_function_2D( mType );
  mf->jacobian( pd, idx, nodes, sample, junk2, junk, junk3, J, err );
  MSQ_ERRZERO(err);
  const MsqVector<3> cross = J.column(0) * J.column(1);

  if (pd.domain_set()) {
    Vector3D norm;
    pd.get_domain_normal_at_sample( pd.get_element_index(this), sample, norm, err );
    MSQ_ERRZERO(err);
  
    const MsqVector<3> N(&norm[0]);
    if (cross % N < 0.0)
      return true;
  }
  
  const double l1 = J.column(0) % J.column(0);
  const double l2 = J.column(1) % J.column(1);
  return cross % cross < DBL_EPSILON*DBL_EPSILON * l1*l2;
}

NodeSet
MsqMeshEntity::all_nodes( MsqError& err ) const
{
  bool mid_edge, mid_face, mid_vol;
  TopologyInfo::higher_order( mType, node_count(), mid_edge, mid_face, mid_vol, err );
  NodeSet result;
  result.set_all_corner_nodes( mType );
  if (mid_edge)
    result.set_all_mid_edge_nodes( mType );
  if (mid_face)
    result.set_all_mid_face_nodes( mType );
  if (mid_vol)
    result.set_all_mid_region_nodes( mType );
  return result;
}

void MsqMeshEntity::check_element_orientation_corners(
  PatchData &pd, 
  int& inverted,
  int& total,
  MsqError &err)
{
  int num_nodes = node_count();
  total = inverted = 0;

  if (node_count() > vertex_count()) {
    MSQ_SETERR(err)("Cannot perform inversion test for higher-order element"
                    " without mapping function.", MsqError::UNSUPPORTED_ELEMENT );
    return;
  }

  const MsqVertex *vertices = pd.get_vertex_array(err);  MSQ_ERRRTN(err);

  const Vector3D d_con(1.0, 1.0, 1.0);

  int i;
  Vector3D coord_vectors[3];
  Vector3D center_vector;
  
  switch(mType) {
    case TRIANGLE:
      
      if (!pd.domain_set())
        return;
      
      pd.get_domain_normal_at_element(this, coord_vectors[2], err); MSQ_ERRRTN(err);
      coord_vectors[2] = coord_vectors[2] / coord_vectors[2].length();// Need unit normal
      center_vector = vertices[vertexIndices[0]];
      coord_vectors[0] = vertices[vertexIndices[1]]-center_vector;
      coord_vectors[1] = vertices[vertexIndices[2]]-center_vector;
      total = 1;
      inverted = (coord_vectors[2]%(coord_vectors[0]*coord_vectors[1] ) <= 0.0);
      break;
    
    case QUADRILATERAL:
      
      if (!pd.domain_set())
        return;

      pd.get_domain_normal_at_element(this, coord_vectors[2], err); MSQ_ERRRTN(err);
      coord_vectors[2] = coord_vectors[2] / coord_vectors[2].length();// Need unit normal

      for (i = 0; i < 4; ++i) {
        center_vector = vertices[vertexIndices[i]];
        coord_vectors[0] = vertices[vertexIndices[(i+1)%4]]-center_vector;
        coord_vectors[1] = vertices[vertexIndices[(i+3)%4]]-center_vector;
        ++total;
        inverted += (coord_vectors[2]%(coord_vectors[0]*coord_vectors[1] ) <= 0.0);
      }
      break;

    case TETRAHEDRON:
      center_vector = vertices[vertexIndices[0]];
      coord_vectors[0] = vertices[vertexIndices[1]]-center_vector;
      coord_vectors[1] = vertices[vertexIndices[2]]-center_vector;
      coord_vectors[2] = vertices[vertexIndices[3]]-center_vector;
      total = 1;
      inverted = ( coord_vectors[0]%(coord_vectors[1]*coord_vectors[2] ) <= 0.0);
      break;

    case POLYGON:
      
      if (!pd.domain_set())
        return;

      pd.get_domain_normal_at_element(this, coord_vectors[2], err); MSQ_ERRRTN(err);
      coord_vectors[2] = coord_vectors[2] / coord_vectors[2].length();// Need unit normal

      for (i = 0; i < num_nodes; ++i) {
        center_vector = vertices[vertexIndices[i]];
        coord_vectors[0] = vertices[vertexIndices[(i+1)%num_nodes]]-center_vector;
        coord_vectors[1] = vertices[vertexIndices[(i+num_nodes-1)%num_nodes]]-center_vector;
        ++total;
        inverted += (coord_vectors[2]%(coord_vectors[0]*coord_vectors[1] ) <= 0.0);
      }
      break;

    default: // generic code for 3D elements
    {
      size_t num_corners = corner_count();
      unsigned num_adj;
      const unsigned* adj_idx;
      for (unsigned j = 0; j < num_corners; ++j)
      {
        adj_idx = TopologyInfo::adjacent_vertices( mType, j, num_adj );
        if (3 != num_adj)
        {
          MSQ_SETERR(err)("Unsupported element type.", MsqError::INVALID_ARG);
          return;
        }
        
        center_vector = vertices[vertexIndices[j]];
        coord_vectors[0] = vertices[vertexIndices[adj_idx[0]]] - center_vector;
        coord_vectors[1] = vertices[vertexIndices[adj_idx[1]]] - center_vector;
        coord_vectors[2] = vertices[vertexIndices[adj_idx[2]]] - center_vector;
        ++total;
        inverted += (coord_vectors[0]%(coord_vectors[1]*coord_vectors[2] ) <= 0.0);
      }
      break;
    }
  } // end switch over element type
}

} // namespace Mesquite
