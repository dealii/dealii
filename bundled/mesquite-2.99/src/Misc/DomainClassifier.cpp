/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2008 Lawrence Livermore National Laboratory.  Under 
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

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */

/** \file DomainClassifier.cpp
 *  \brief 
 *  \author Jason Kraftcheck
 */

#include "DomainClassifier.hpp"
#include "MsqError.hpp"
#include "TopologyInfo.hpp"
#include "MsqVertex.hpp"
#include <map>
#include <set>
#include <algorithm>
#include <iterator>

namespace MESQUITE_NS {

  /**\brief Helper function for 'find_skin'
   *
   * Check the vertices on element sides to see if the
   * elements are adjacent.
   *\param num_vtx  Number of vertices in each element side
   *\param idx_set_1 Indices into the first element's vertex list
   *                 at which the side vertices are located.
   *\param vtx_set_1 Vertex of first element.
   *\param idx_set_2 Indices into the second element's vertex list
   *                 at which the side vertices are located.
   *\param vtx_set_2 Vertex of second element.
   */
static bool compare_sides( unsigned num_vtx,
                           const unsigned* idx_set_1,
                           const Mesh::VertexHandle* vtx_set_1,
                           const unsigned* idx_set_2,
                           const Mesh::VertexHandle* vtx_set_2 );

  /**\brief Skin mesh
   *
   * Find vertices and lower-dimension elements on the boundary
   * of a mesh.  
   */
 static void find_skin( Mesh* mesh,
                        std::vector<Mesh::VertexHandle>& skin_verts,
                        std::vector<Mesh::ElementHandle>& skin_elems,
                        MsqError& err );

  /**\brief Group MeshDomain objectects by topological dimension.
   *
   * Given an array of MeshDomain pointers and and array of dimensions,
   * where each dimension value is the topological dimension of the
   * corresponding domain, construct a DomainSet vector containing
   * the domains, grouped from lowest dimension to highest.
   *\param domains  Array of domains
   *\param dims     Topological dimension of each domain
   *\param count    Length of 'domains' and 'dims' arrays
   *\param results  Output
   *\param dim_indices Output: index in output at which the first
   *                domain of the corresponding topological dimension
   *                occurs.
   */
static void dimension_sort_domains( MeshDomain** domains,
                                    const int* dims,
                                    unsigned count,
                                    std::vector<DomainClassifier::DomainSet>& results,
                                    int dim_indices[4],
                                    MsqError& err );

  /**\brief Classify mesh vertices using vertex position
   *
   * Classify vertices by distance from each geometric domain.
   * A vertex is assigned to a domain if it lies within the domain.
   * If a vertex lies in multiple domains, it is assigned to the
   * domain with the lowest topological dimension.
   *\param mesh  The Mesquite::Mesh instance
   *\param dim_sorted_domains  The list of MeshDomains, sorted from lowest 
   *             to highest dimension.
   *\param num_domain  Length of 'dim_sorted_domains' array.
   *\param vertices    Vertices to classify
   *\param epsilon     Maximum distance a vertex may deviate from its domain.
   */
static void geom_classify_vertices( Mesh* mesh,
                                    DomainClassifier::DomainSet* dim_sorted_domains,
                                    unsigned num_domain,
                                    const std::vector<Mesh::VertexHandle>& vertices,
                                    double epsilon,
                                    MsqError& err );

  /**\brief Classify mesh elements using vertex classification
   *
   * Classify elements using vertex classification, according to the
   * following rules, listed highest-precedence first.
   * Case 1: one or more vertices has no classification -> element has no classification
   * Case 2: all vertices classified to domain with dimension < 2 -> unknown
   * Case 3: one or more vertices has dimension-2 classification:
   *  Case 3a: all have same dimension-2 domain -> element assigned to same domain
   *  Case 3b: different dimension-2 domains -> element has no classification
   *
   *\param mesh  The Mesquite::Mesh instance
   *\param dim_sorted_domains  The list of MeshDomains, sorted from lowest 
   *             to highest dimension.
   *\param num_domain  Length of 'dim_sorted_domains' array.
   *\param elems    Elements to classify
   *\param uknown   Elements for which all vertices are classified to 
   *                domains with dimension less than two.
   */
static void vert_classify_elements( Mesh* mesh,
                                    DomainClassifier::DomainSet* dim_sorted_domains,
                                    unsigned num_domain,
                                    int dim_indices[4],
                                    const std::vector<Mesh::ElementHandle>& elems,
                                    std::vector<Mesh::ElementHandle>& unknown_elems,
                                    MsqError& err );

  /**\brief Classify mesh elements geometrically. */
static void geom_classify_elements( Mesh* mesh,
                                    DomainClassifier::DomainSet* dim_sorted_domains,
                                    unsigned num_domain,
                                    int dim_indices[4],
                                    const std::vector<Mesh::ElementHandle>& elems,
                                    std::vector<Mesh::ElementHandle>& unknown_elems,
                                    double epsilon,
                                    MsqError& err );
                                    
                                    
                                    
static bool compare_sides( unsigned num_vtx,
                           const unsigned* idx_set_1,
                           const Mesh::VertexHandle* vtx_set_1,
                           const unsigned* idx_set_2,
                           const Mesh::VertexHandle* vtx_set_2 )
{
  unsigned i;
  for (i = 0;; ++i) {
    if (i >= num_vtx)
      return false;
    if (vtx_set_1[idx_set_1[0]] == vtx_set_2[idx_set_2[i]])
      break;
  }
  
  if (num_vtx == 1)
    return true;
  
  if (vtx_set_1[idx_set_1[1]] == vtx_set_2[idx_set_2[(i+1)%num_vtx]]) {
    for (unsigned j = 2; j < num_vtx; ++j)
      if (vtx_set_1[idx_set_1[j]] != vtx_set_2[idx_set_2[(i+j)%num_vtx]])
        return false;
    return true;
  }
  else if (vtx_set_1[idx_set_1[1]] == vtx_set_2[idx_set_2[(i+num_vtx-1)%num_vtx]]) {
    for (unsigned j = 2; j < num_vtx; ++j)
      if (vtx_set_1[idx_set_1[j]] != vtx_set_2[idx_set_2[(i+num_vtx-j)%num_vtx]])
        return false;
    return true;
  }
  else {
    return false;
  }
}

static void find_skin( Mesh* mesh,
                       std::vector<Mesh::VertexHandle>& skin_verts,
                       std::vector<Mesh::ElementHandle>& skin_elems,
                       MsqError& err )
{
  std::vector<Mesh::ElementHandle> elements, adj_elem, side_elem;
  std::vector<Mesh::VertexHandle> vertices, adj_vtx;
  std::vector<size_t> junk;
  mesh->get_all_elements( elements, err ); MSQ_ERRRTN(err);
  if (elements.empty())
    return;
  const unsigned elem_idx[] = { 0, 1, 2, 3 };
  
  std::vector<unsigned> extra_side_vtx;
  for (size_t e = 0; e < elements.size(); ++e) {
    Mesh::ElementHandle elem = elements[e];
    vertices.clear();
    mesh->elements_get_attached_vertices( &elem, 1, vertices, junk, err ); MSQ_ERRRTN(err);

    EntityTopology type;
    mesh->elements_get_topologies( &elem, &type, 1, err ); MSQ_ERRRTN(err);
    const unsigned dim = TopologyInfo::dimension( type );
    const unsigned nsides = TopologyInfo::sides( type );
    for (unsigned s = 0; s < nsides; ++s) {
      unsigned num_side_vtx;
      const unsigned* side_vtx = TopologyInfo::side_vertices( type, dim-1, s, num_side_vtx );
      adj_elem.clear();
      side_elem.clear();
      mesh->vertices_get_attached_elements( &vertices[*side_vtx], 1, adj_elem, junk, err ); MSQ_ERRRTN(err);
      
      bool found_adj = false;
      for (size_t a = 0; a < adj_elem.size(); ++a) {
        if (adj_elem[a] == elem)
          continue;
        EntityTopology adj_type;
        mesh->elements_get_topologies( &adj_elem[a], &adj_type, 1, err ); MSQ_ERRRTN(err);
        adj_vtx.clear();
        mesh->elements_get_attached_vertices( &adj_elem[a], 1, adj_vtx, junk, err ); MSQ_ERRRTN(err);
        if (TopologyInfo::dimension(adj_type) == dim) {
          unsigned nsides2= TopologyInfo::sides( adj_type );
          for (unsigned s2 = 0; s2 < nsides2; ++s2) {
            unsigned num_side_vtx2;
            const unsigned* side_vtx2 = TopologyInfo::side_vertices( adj_type, dim-1, s2, num_side_vtx2 );
            if (num_side_vtx2 == num_side_vtx &&
                compare_sides( num_side_vtx, side_vtx, arrptr(vertices), side_vtx2, arrptr(adj_vtx) )) {
              found_adj = true;
              break;
            }
          }
        }
        else {
          if (TopologyInfo::corners(adj_type) == num_side_vtx)
            if (compare_sides( num_side_vtx, side_vtx, arrptr(vertices), elem_idx, arrptr(adj_vtx) ))
              side_elem.push_back( adj_elem[a] );
        }
      }
      
      if (!found_adj) {
        int ho = TopologyInfo::higher_order( type, vertices.size(), err );  MSQ_ERRRTN(err);
        // mask mid-element node (because that cannot be on the skin)
        ho &= ~(1<<dim);
        // if element has other higher-order nodes, we need to find the subset on the skin side
        if (ho) {
          extra_side_vtx.clear();
          extra_side_vtx.insert( extra_side_vtx.end(), side_vtx, side_vtx+num_side_vtx );
            // get mid-side node
          if (ho & (1<<(dim-1))) {
            int idx = TopologyInfo::higher_order_from_side( type, vertices.size(), dim-1, s, err ); MSQ_ERRRTN(err);
            extra_side_vtx.push_back( idx );
          }
          if (dim == 3 && (ho & 2)) { // if volume element, need mid-edge nodes too
            const unsigned nedges = TopologyInfo::edges( type );
            for (unsigned e = 0; e < nedges; ++e) {
              unsigned junk;
              const unsigned* edge_vtx = TopologyInfo::side_vertices( type, 1, e, junk );
              assert(edge_vtx && 2 == junk);
              unsigned idx = std::find( side_vtx, side_vtx+num_side_vtx, edge_vtx[0] ) - side_vtx;
              if (idx < num_side_vtx) {
                if ((edge_vtx[1] == side_vtx[(idx+1)%num_side_vtx]) ||
                    (edge_vtx[1] == side_vtx[(idx+num_side_vtx-1)%num_side_vtx])) {
                  idx = TopologyInfo::higher_order_from_side( type, vertices.size(), 1, e, err ); MSQ_ERRRTN(err);
                  extra_side_vtx.push_back( idx );
                }
              }
            }
          }
          
          num_side_vtx = extra_side_vtx.size();
          side_vtx = arrptr(extra_side_vtx);
        }
      
        for (unsigned v = 0; v < num_side_vtx; ++v)
          skin_verts.push_back( vertices[side_vtx[v]] );
        for (unsigned j = 0; j < side_elem.size(); ++j)
          skin_elems.push_back( side_elem[j] );
      }
    }
  }
  
  std::sort( skin_verts.begin(), skin_verts.end() );
  skin_verts.erase( std::unique( skin_verts.begin(), skin_verts.end() ), skin_verts.end() );
  std::sort( skin_elems.begin(), skin_elems.end() );
  skin_elems.erase( std::unique( skin_elems.begin(), skin_elems.end() ), skin_elems.end() );
}

static void dimension_sort_domains( MeshDomain** domains,
                                    const int* dims,
                                    unsigned count,
                                    std::vector<DomainClassifier::DomainSet>& results,
                                    int dim_indices[4],
                                    MsqError& err )
{
  results.clear();
  results.resize(count);
  unsigned idx = 0;
  int dim;
  for (dim = 0; idx < count; ++dim) {
    if (dim > 3) {
      MSQ_SETERR(err)("Invalid domain dimension", MsqError::INVALID_ARG);
      return;
    }
    dim_indices[dim] = idx;
    for (unsigned j = 0; j < count; ++j) {
      if (dims[j] == dim) {
        results[idx].domain = domains[j];
        ++idx;
      }
    }
  }
  for (; dim <= 3; ++dim)
    dim_indices[dim] = idx;
}

static void geom_classify_vertices( Mesh* mesh,
                                    DomainClassifier::DomainSet* dim_sorted_domains,
                                    unsigned num_domain,
                                    const std::vector<Mesh::VertexHandle>& vertices,
                                    double epsilon,
                                    MsqError& err )
{
  Vector3D pt;
  MsqVertex coord;
  unsigned i;
  unsigned short dim;
  const double epssqr = epsilon*epsilon;
  std::vector<Mesh::VertexHandle>::const_iterator iter;
  for (iter = vertices.begin(); iter != vertices.end(); ++iter) {
    mesh->vertices_get_coordinates( &*iter, &coord, 1, err ); MSQ_ERRRTN(err);
    for (i = 0; i < num_domain; ++i) {
        // Call DoF function first just to see if it fails
        // (e.g. if domain knows which vertices it is responsible for,
        //  it should fail if this vertex isn't one.)
      dim_sorted_domains[i].domain->domain_DoF( &*iter, &dim, 1 , err );
      if (err || dim > 2) {
        err.clear();
        continue;
      }
      pt = coord;
      dim_sorted_domains[i].domain->snap_to( *iter, pt );
      if ((pt - coord).length_squared() <= epssqr) {
        dim_sorted_domains[i].vertices.push_back( *iter );
        break;
      }
    }
  }
}

static void vert_classify_elements( Mesh* mesh,
                                    DomainClassifier::DomainSet* dim_sorted_domains,
                                    unsigned num_domain,
                                    int dim_indices[4],
                                    const std::vector<Mesh::ElementHandle>& elems,
                                    std::vector<Mesh::ElementHandle>& unknown_elems,
                                    MsqError& err )
{
    // sort vertex lists for faster search
  for (unsigned i = 0; i < num_domain; ++i)
    std::sort( dim_sorted_domains[i].vertices.begin(), dim_sorted_domains[i].vertices.end() );
  
  std::vector<Mesh::ElementHandle>::const_iterator iter;
  std::vector<Mesh::VertexHandle> verts;
  std::vector<size_t> junk;
  for (iter = elems.begin(); iter != elems.end(); ++iter) {
    int dom = -1;
    verts.clear();
    mesh->elements_get_attached_vertices( &*iter, 1, verts, junk, err ); MSQ_ERRRTN(err);
    for (unsigned v = 0; v < verts.size(); ++v) {
      int i;
      for (i = 0; i < dim_indices[3]; ++i) {
        DomainClassifier::DomainSet* d = dim_sorted_domains + i;
        if (std::binary_search(d->vertices.begin(), d->vertices.end(), verts[v]))
          break;
      }
        // if any vertex in element has no domain, then element has no domain
      if (i < 0) {
        dom = -2; 
        break;
      }
        // if vertex is in curve or point, ignore it
      else if (i < dim_indices[2])
        continue;
        // else if we already have a vertex on a different surface,
        // element must be in volume (no domain)
      else if (dom >= 0 && dom != i) {
        dom = -2; 
        break;
      }
      else {
        dom = i;
      }
    }
    
    if (dom >= 0)
      dim_sorted_domains[dom].elements.push_back( *iter );
    else if (dom == -1) // all verts on curves or points
      unknown_elems.push_back( *iter );
  }
}

static void geom_classify_elements( Mesh* mesh,
                                    DomainClassifier::DomainSet* dim_sorted_domains,
                                    unsigned num_domain,
                                    int dim_indices[4],
                                    const std::vector<Mesh::ElementHandle>& elems,
                                    std::vector<Mesh::ElementHandle>& unknown_elems,
                                    double epsilon,
                                    MsqError& err )
{
  if (elems.empty())
    return;

  const double epssqr = epsilon*epsilon;
  Vector3D pt;
  MsqVertex coord;
  std::vector<Mesh::ElementHandle>::const_iterator iter;
  std::vector<Mesh::VertexHandle> verts;
  std::vector<MsqVertex> coords;
  std::vector<size_t> junk;
  for (iter = elems.begin(); iter != elems.end(); ++iter) {
    verts.clear();
    mesh->elements_get_attached_vertices( &*iter, 1, verts, junk, err ); MSQ_ERRRTN(err);
    coords.resize(verts.size());
    mesh->vertices_get_coordinates( arrptr(verts), arrptr(coords), verts.size(), err ); MSQ_ERRRTN(err);

      // get all 2D domains that contain first vertex
    std::vector<int> doms;
    for (int i = dim_indices[2]; i < dim_indices[3]; ++i) {
      pt = coords[0];
      dim_sorted_domains[i].domain->snap_to( verts[0], pt );
      if ((pt - coords[0]).length_squared() <= epssqr)
        doms.push_back( i );
    }
    
      // of the 2D domains containing the first vertex,
      // remove any that do not contain all other corners of the element.
    for (unsigned i = 1; i < verts.size(); ++i) {
      for (unsigned j = 0; j < doms.size(); ) {
        pt = coords[i];
        dim_sorted_domains[doms[j]].domain->snap_to( verts[j], pt );
        if ((pt - coords[0]).length_squared() <= epssqr)
          ++j;
        else
          doms.erase( doms.begin()+j );
      }
    }
    
    if (doms.size() == 1) 
      dim_sorted_domains[doms.front()].elements.push_back( *iter );
    else
      unknown_elems.push_back( *iter );
  }
}

void DomainClassifier::classify_geometrically( DomainClassifier& result,
                                        Mesh* mesh,
                                        double tolerance,
                                        MeshDomain** domain_array,
                                        const int* dimension_array,
                                        unsigned array_length,
                                        MsqError& err )
{
  if (!array_length)
    return;

  std::vector<DomainSet> domains;
  int dim_indices[4];
  dimension_sort_domains( domain_array, dimension_array, array_length, domains, dim_indices, err );
  MSQ_ERRRTN(err);
  
    // classify vertices by position
  std::vector<Mesh::VertexHandle> vertices;
  mesh->get_all_vertices( vertices, err ); MSQ_ERRRTN(err);
  geom_classify_vertices( mesh, arrptr(domains), dim_indices[3], vertices, tolerance, err );
  MSQ_ERRRTN(err);
  
    // get elements and types
  std::vector<Mesh::ElementHandle> elems;
  mesh->get_all_elements( elems, err ); MSQ_ERRRTN(err);
  if (elems.empty())
    return;
  std::vector<EntityTopology> types( elems.size() );
  mesh->elements_get_topologies( arrptr(elems), arrptr(types), elems.size(), err );
  MSQ_ERRRTN(err);
  
    // get rid of elements w/ dimension other than 2
  size_t w = 0;
  for (size_t r = 0; r < elems.size(); ++r) {
    if (TopologyInfo::dimension(types[r]) == 2) {
      elems[w] = elems[r];
      ++w;
    }
  }
  elems.resize( w );
 
    // classify elements using vertex classification
  std::vector<Mesh::ElementHandle> unknown;
  vert_classify_elements( mesh, arrptr(domains), domains.size(), dim_indices, elems, unknown, err );
  MSQ_ERRRTN(err);
  
    // classify unknown elements geometrically
  elems.swap( unknown );
  unknown.clear();
  geom_classify_elements( mesh, arrptr(domains), domains.size(), dim_indices, elems, unknown, tolerance, err );
  
  classify_by_handle( result, mesh, arrptr(domains), domains.size(), err );
}
  

void DomainClassifier::classify_skin_geometrically( DomainClassifier& result,
                                        Mesh* mesh,
                                        double tolerance,
                                        MeshDomain** domain_array,
                                        const int* dimension_array,
                                        unsigned array_length,
                                        MsqError& err )
{
  if (!array_length)
    return;

  std::vector<DomainSet> domains;
  int dim_indices[4];
  dimension_sort_domains( domain_array, dimension_array, array_length, domains, dim_indices, err );
  MSQ_ERRRTN(err);
  
  std::vector<Mesh::VertexHandle> vertices;
  std::vector<Mesh::ElementHandle> elements;
  find_skin( mesh, vertices, elements, err ); MSQ_ERRRTN(err);
  
    // classify vertices by position
  geom_classify_vertices( mesh, arrptr(domains), dim_indices[3], vertices, tolerance, err );
  MSQ_ERRRTN(err);
 
    // classify elements using vertex classification
  std::vector<Mesh::ElementHandle> unknown;
  vert_classify_elements( mesh, arrptr(domains), domains.size(), dim_indices, elements, unknown, err );
  MSQ_ERRRTN(err);
  
    // classify unknown elements geometrically
  elements.swap( unknown );
  unknown.clear();
  geom_classify_elements( mesh, arrptr(domains), domains.size(), dim_indices, elements, unknown, tolerance, err );

  if (!unknown.empty()) {
    MSQ_SETERR(err)("Count not classify all skin elements", MsqError::INVALID_MESH );
    return;
  }
  
  classify_by_handle( result, mesh, arrptr(domains), domains.size(), err );
}
  


void DomainClassifier::classify_by_tag( DomainClassifier& result,
                                        Mesh* mesh,
                                        const char* tag_name,
                                        MeshDomain** domain_array,
                                        const int* id_array,
                                        unsigned array_length,
                                        MsqError& err )
{
  TagHandle tag = mesh->tag_get( tag_name, err ); MSQ_ERRRTN(err);
  
  std::string name2;
  Mesh::TagType type;
  unsigned size;
  mesh->tag_properties( tag, name2, type, size, err ); MSQ_ERRRTN(err);
  if (type != Mesh::INT || size != 1) {
    MSQ_SETERR(err)(MsqError::TAG_NOT_FOUND,"Tag does not contain single integer value: %s\n", tag_name );
    return;
  }
  
  std::vector<DomainSet> sets(array_length);
  std::map<int,DomainSet*> idmap;
  for (unsigned i = 0; i < array_length; ++i) {
    sets[i].domain = domain_array[i];
    idmap[id_array[i]] = &sets[i];
  }
  
  std::vector<Mesh::VertexHandle> vertices;
  std::vector<Mesh::ElementHandle> elements;
  mesh->get_all_vertices( vertices, err ); MSQ_ERRRTN(err);
  mesh->get_all_elements( elements, err ); MSQ_ERRRTN(err);
  if (vertices.empty())
    return;
  
  std::map<int,DomainSet*>::const_iterator iter;
  std::vector<int> ids( vertices.size() );
  mesh->tag_get_vertex_data( tag, vertices.size(), arrptr(vertices), arrptr(ids), err );  MSQ_ERRRTN(err);
  for (size_t j = 0; j < vertices.size(); ++j) {
    iter = idmap.find( ids[j] );
    if (iter != idmap.end())
      iter->second->vertices.push_back( vertices[j] );
  }
  
  ids.clear();
  ids.resize( elements.size() );
  if (!elements.empty()) {
    mesh->tag_get_element_data( tag, elements.size(), arrptr(elements), arrptr(ids), err ); 
    if (err) {
      err.clear();
    }
    else {
      for (size_t k = 0; k < elements.size(); ++k) {
        iter = idmap.find( ids[k] );
        if (iter != idmap.end())
          iter->second->elements.push_back( elements[k] );
      }
    }
  }
  
  if (!sets.empty())
    classify_by_handle( result, mesh, arrptr(sets), sets.size(), err );
}
                              

void DomainClassifier::classify_by_handle( DomainClassifier& result,
                                           Mesh* mesh,
                                           DomainSet* sets,
                                           unsigned array_length,
                                           MsqError& err )
{
  result.clear();

    // get all vertices and elements
  std::vector<Mesh::VertexHandle> vertices;
  std::vector<Mesh::ElementHandle> elements;
  mesh->get_all_vertices( vertices, err ); MSQ_ERRRTN(err);
  mesh->get_all_elements( elements, err ); MSQ_ERRRTN(err);

    // sort all arrays so we can merge
  std::sort( vertices.begin(), vertices.end() );
  std::sort( elements.begin(), elements.end() );
  for (unsigned i = 0; i < array_length; ++i) {
    std::sort( sets[i].vertices.begin(), sets[i].vertices.end() );
    std::sort( sets[i].elements.begin(), sets[i].elements.end() );
  }
  
    // build vertex block list
  std::vector<size_t> indices(array_length, 0);
  size_t idx = 0;
  while (idx < vertices.size()) {
    DomainBlock block;
    block.firstHandle = vertices[idx];
    block.lastHandle = vertices[idx];
    
      // find domain
    unsigned dom = 0;
    for (dom = 0; dom < array_length; ++dom) {
      const size_t i = indices[dom];
      if (i < sets[dom].vertices.size() && sets[dom].vertices[i] == vertices[idx]) {
        ++indices[dom];
        break;
      }
    }
    ++idx;
    if (dom == indices.size()) // no domain
      continue;
      
    block.domain = sets[dom].domain;
    while (indices[dom] < sets[dom].vertices.size() &&
           sets[dom].vertices[indices[dom]] == vertices[idx]) {
      block.lastHandle = vertices[idx];
      ++idx;
      ++indices[dom];
    }
    result.vertexList.push_back( block );
  }
  
    // build vertex block list
  indices.clear();
  indices.resize( array_length, 0 );
  idx = 0;
  while (idx < elements.size()) {
    DomainBlock block;
    block.firstHandle = elements[idx];
    block.lastHandle = elements[idx];
    
      // find domain
    unsigned dom = 0;
    for (dom = 0; dom < array_length; ++dom) {
      const size_t i = indices[dom];
      if (i < sets[dom].elements.size() && sets[dom].elements[i] == elements[idx]) {
        ++indices[dom];
        break;
      }
    }
    ++idx;
    if (dom == indices.size()) // no domain
      continue;
      
    block.domain = sets[dom].domain;
    while (indices[dom] < sets[dom].elements.size() &&
           sets[dom].elements[indices[dom]] == elements[idx]) {
      block.lastHandle = elements[idx];
      ++idx;
      ++indices[dom];
    }
    result.elementList.push_back( block );
  }
/*  
  printf("\n");
  for (unsigned i = 0; i < result.vertexList.size(); ++i) 
    printf("v %2lu-%2lu @ %p\n", (unsigned long)result.vertexList[i].firstHandle,
                               (unsigned long)result.vertexList[i].lastHandle,
                               result.vertexList[i].domain );
  for (unsigned i = 0; i < result.elementList.size(); ++i) 
    printf("e %2lu-%2lu @ %p\n", (unsigned long)result.elementList[i].firstHandle,
                               (unsigned long)result.elementList[i].lastHandle,
                               result.elementList[i].domain );
  printf("\n");
*/
}

static bool next_vertex( Mesh* mesh, 
                         Mesh::VertexHandle& vtx,
                         std::set<Mesh::VertexHandle>& unseen,
                         MsqError& err )
{
  std::vector<Mesh::ElementHandle> vtx_elems;
  std::vector<size_t> junk;
  
  mesh->vertices_get_attached_elements( &vtx, 1, vtx_elems, junk, err );
  MSQ_ERRZERO(err);

  std::vector<EntityTopology> elem_types( vtx_elems.size() );
  if (!vtx_elems.empty()) {
    mesh->elements_get_topologies( arrptr(vtx_elems), arrptr(elem_types), vtx_elems.size(), err );
    MSQ_ERRZERO(err);
  }

  std::vector<Mesh::VertexHandle> corners;
  for (size_t j = 0; j < vtx_elems.size(); ++j) {
    corners.clear();
    mesh->elements_get_attached_vertices( &vtx_elems[j], 1, corners, junk, err );
    MSQ_ERRZERO(err);

    unsigned nedges = TopologyInfo::edges( elem_types[j] );
    unsigned vidx = std::find( corners.begin(), corners.end(), vtx ) - corners.begin();
    

    // Check mid-edge nodes first, if present
    int ho = TopologyInfo::higher_order( elem_types[j], corners.size(), err );
    MSQ_ERRZERO(err);
      // If element has mid-edge nodes *and* current vertex is a corner vertex
    if ((ho & 2) && (vidx < TopologyInfo::corners(elem_types[j]))) {
      for (unsigned e = 0; e < nedges; ++e) {
        const unsigned* edge_verts = TopologyInfo::edge_vertices( elem_types[j], e, err );
        MSQ_ERRZERO(err);
        if (edge_verts[0] == vidx || edge_verts[1] == vidx) {
          int idx = TopologyInfo::higher_order_from_side( elem_types[j], corners.size(), 1, e, err );
          MSQ_ERRZERO(err);
          std::set<Mesh::VertexHandle>::iterator f = unseen.find( corners[idx] );
          if (f != unseen.end()) {
            vtx = *f;
            unseen.erase(f);
            return true;
          }
        }
      }
    }

      // if current vertx is a mid-edge node
    else if (ho & 2) {
      unsigned d, e;
      TopologyInfo::side_from_higher_order( elem_types[j], corners.size(), vidx, d, e, err );
      MSQ_ERRZERO(err);
      if (d != 1)
        continue;
      const unsigned* edge_verts =  TopologyInfo::edge_vertices( elem_types[j], e, err );
      MSQ_ERRZERO(err);
      for (int v = 0; v < 2; ++v) {
        std::set<Mesh::VertexHandle>::iterator f = unseen.find( corners[edge_verts[v]] );
        if (f != unseen.end()) {
          vtx = *f;
          unseen.erase(f);
          return true;
        }
      }
    }
    
    else {
      for (unsigned e = 0; e < nedges; ++e) {
        const unsigned* edge_verts = TopologyInfo::edge_vertices( elem_types[j], e, err );
        MSQ_ERRZERO(err);
        int idx;
        if (edge_verts[0] == vidx)
          idx = edge_verts[1];
        else if (edge_verts[1] == vidx) 
          idx = edge_verts[0];
        else
          continue;
          
        std::set<Mesh::VertexHandle>::iterator f = unseen.find( corners[idx] );
        if (f != unseen.end()) {
          vtx = *f;
          unseen.erase(f);
          return true;
        }
      }
    }
  }
  
  return false;
} 
  

void DomainClassifier::test_valid_classification( Mesh* mesh,
                                                  MsqError& err )
{
  size_t i;

  // Get all mesh entities
  std::vector<Mesh::VertexHandle> vertices;
  std::vector<Mesh::ElementHandle> elements;
  mesh->get_all_vertices( vertices, err ); MSQ_ERRRTN(err);
  mesh->get_all_elements( elements, err ); MSQ_ERRRTN(err);
  std::sort( vertices.begin(), vertices.end() );
  std::sort( elements.begin(), elements.end() );
  
  // Get contents of each domain.
  std::map<MeshDomain*,int>::iterator iter;
  std::map<MeshDomain*,int> idxmap;
  for (i = 0; i < vertexList.size(); ++i)
    idxmap[vertexList[i].domain] = 0;
  for (i = 0; i < elementList.size(); ++i)
    idxmap[elementList[i].domain] = 0;
  std::vector<DomainSet> domains(idxmap.size());
  int idx = 0;
  for (iter = idxmap.begin(); iter != idxmap.end(); ++iter) {
    iter->second = idx;
    domains[idx].domain = iter->first;
    ++idx;
  }
  for (i = 0; i < vertexList.size(); ++i) {
    std::vector<Mesh::VertexHandle>::const_iterator s, e;
    s = std::lower_bound( vertices.begin(), vertices.end(), vertexList[i].firstHandle );
    e = std::upper_bound( vertices.begin(), vertices.end(), vertexList[i].lastHandle );
    DomainSet& set = domains[idxmap[vertexList[i].domain]];
    std::copy( s, e, std::back_inserter(set.vertices) );
  }
  for (i = 0; i < elementList.size(); ++i) {
    std::vector<Mesh::ElementHandle>::const_iterator s, e;
    s = std::lower_bound( elements.begin(), elements.end(), elementList[i].firstHandle );
    e = std::upper_bound( elements.begin(), elements.end(), elementList[i].lastHandle );
    DomainSet& set = domains[idxmap[elementList[i].domain]];
    std::copy( s, e, std::back_inserter(set.elements) );
  }

  // guess geometric dimension for each domain
  std::map<DomainSet*,int> dimmap;
  std::vector<unsigned short> dof;
  for (i = 0; i < domains.size(); ++i) {
    if (!domains[i].elements.empty())
      dimmap[&domains[i]] = 2;
    else {
      dof.resize( domains[i].vertices.size() );
      domains[i].domain->domain_DoF( &(domains[i].vertices[0]), arrptr(dof), dof.size(), err );
      MSQ_ERRRTN(err);
      unsigned short dim = *std::max_element( dof.begin(), dof.end() );
      dimmap[&domains[i]] = dim;
    }
  }
  
  // group domains by dimension
  std::vector<DomainSet*> points, curves, surfaces;
  for (std::map<DomainSet*,int>::iterator it = dimmap.begin(); it != dimmap.end(); ++it) {
    switch (it->second) {
      case 0:
        points.push_back( it->first );
        break;
      case 1:
        curves.push_back( it->first );
        break;
      case 2:
        surfaces.push_back( it->first );
        break;
      default:
        MSQ_SETERR(err)(MsqError::INVALID_STATE, "Invalid domain dimension: %d\n", it->second);
        break;
    }
  }
  
  // check that each point-domain has a single vertex
  for (i = 0; i < points.size(); ++i) {
    if (points[i]->vertices.size() != 1 ||
       !points[i]->elements.empty()) {
      MSQ_SETERR(err)("Point domain mesh not a single vertex.", MsqError::INVALID_STATE);
      return;
    }
  }
  
  // check that each curve domain has a chain of connected edges
  std::set<Mesh::VertexHandle> unseen;
  for (i = 0; i < curves.size(); ++i) {
    if (!curves[i]->elements.empty()) {
      MSQ_SETERR(err)("Elements associated with 1D domain.",MsqError::INVALID_STATE);
      return;
    }
  
    unseen.clear();
    std::copy( curves[i]->vertices.begin(),
                   curves[i]->vertices.end(),
                   std::inserter(unseen,unseen.begin()) );
    
    const Mesh::VertexHandle first_vtx = *unseen.begin();
    unseen.erase(unseen.begin());
    Mesh::VertexHandle vtx = first_vtx; 
      // find chain of vertices
    while (next_vertex(mesh, vtx, unseen, err))
      MSQ_ERRRTN(err);
      // search again from the starting vertex because
      // it probably wasn't the first one on the curve
    vtx = first_vtx;
    while (next_vertex(mesh, vtx, unseen, err))
      MSQ_ERRRTN(err);
      // were all vertices in a chain?
    if (!unseen.empty()) {
      MSQ_SETERR(err)("Curve domain contains vertices not in a simply connected chain.", MsqError::INVALID_STATE );
      return;
    }
  }
  
  
  std::set<Mesh::VertexHandle> seen;
  std::set<Mesh::ElementHandle> remaining;
  std::vector<Mesh::VertexHandle> verts, verts2;
  std::vector<Mesh::ElementHandle> stack, vert_elems;
  std::vector<EntityTopology> types;
  std::vector<size_t> junk;
    // if surface contains elements...
  for (i = 0; i < surfaces.size(); ++i) {
    if (surfaces[i]->elements.empty())
      continue;
      
      // Check that any vertices on surface are contained in an
      // element on the surface.
    verts.clear();
    mesh->elements_get_attached_vertices( &(surfaces[i]->elements[0]),
                                          surfaces[i]->elements.size(),
                                          verts, junk, err ); MSQ_ERRRTN(err);
    seen.clear();
    std::copy( verts.begin(), verts.end(), std::inserter(seen,seen.begin()) );

    std::vector<Mesh::VertexHandle>::const_iterator v;
    for (v = surfaces[i]->vertices.begin(); v != surfaces[i]->vertices.end(); ++v) {
      std::set<Mesh::VertexHandle>::iterator j = seen.find(*v);
      if (j == seen.end()) {
        MSQ_SETERR(err)("Vertex on surface domain not in any element.",
                        MsqError::INVALID_STATE);
        return;
      }
    }
      
      // check that elements form 2D patch
    stack.clear();
    remaining.clear();
    std::copy( surfaces[i]->elements.begin(),
                   surfaces[i]->elements.end(),
                   std::inserter(remaining,remaining.begin()) );
    stack.push_back( *remaining.begin() );
    remaining.erase(remaining.begin());
    while (!stack.empty()) {
      Mesh::ElementHandle elem = stack.back();
      stack.pop_back();
      verts.clear();
      mesh->elements_get_attached_vertices( &elem, 1, verts, junk, err ); MSQ_ERRRTN(err);
      // for each edge
      for (size_t j = 0; j < verts.size(); ++j) {
        Mesh::VertexHandle v1 = verts[j], v2 = verts[(j+1)%verts.size()];
        vert_elems.clear();
        mesh->vertices_get_attached_elements( &v1, 1, vert_elems, junk, 
                                              err ); MSQ_ERRRTN(err);
        types.resize( vert_elems.size() );
        if (!vert_elems.empty()) {
          mesh->elements_get_topologies( arrptr(vert_elems), arrptr(types), 
                                 vert_elems.size(), err ); MSQ_ERRRTN(err);
        }
        while (!vert_elems.empty()) {
          Mesh::ElementHandle e2 = vert_elems.back();
          EntityTopology type = types.back();
          vert_elems.pop_back();
          types.pop_back();
          if (TopologyInfo::dimension(type) != 2)
            continue;
          verts2.clear();
          mesh->elements_get_attached_vertices( &e2, 1, verts2, junk, err ); MSQ_ERRRTN(err);
          size_t idx = std::find(verts2.begin(), verts2.end(), v1 ) - verts2.begin();
          if (verts2[(idx+1)%verts2.size()] != v2 &&
              verts2[(idx+verts2.size()-1)%verts2.size()] != v2)
            continue;
          std::set<Mesh::ElementHandle>::iterator r = remaining.find(e2);
          if (r == remaining.end())
            continue;

          stack.push_back(*r);
          remaining.erase(r);
        }
      }
    }

    if (!remaining.empty()) {
      MSQ_SETERR(err)("Surface mesh not a single, simply-connected patch", 
                      MsqError::INVALID_STATE);
      return;
    }
  }
    
      // check that sides of volume elements that are on surface
      // form simply connected patch
  for (i = 0; i < surfaces.size(); ++i) {
      // build list of sides known to be on surface
    std::sort( surfaces[i]->vertices.begin(), surfaces[i]->vertices.end() );
    std::set< std::pair<Mesh::ElementHandle,int> > sides;
    for (size_t j= 0; j < surfaces[i]->vertices.size(); ++j) {
      Mesh::VertexHandle v = surfaces[i]->vertices[j];
      vert_elems.clear();
      mesh->vertices_get_attached_elements( &v, 1, vert_elems, junk, 
                                            err ); MSQ_ERRRTN(err);
      types.resize( vert_elems.size() );
      if (!vert_elems.empty()) {
        mesh->elements_get_topologies( arrptr(vert_elems), arrptr(types), 
                               vert_elems.size(), err ); MSQ_ERRRTN(err);
      }
      while (!vert_elems.empty()) {
        Mesh::ElementHandle e = vert_elems.back();
        EntityTopology type = types.back();
        vert_elems.pop_back();
        types.pop_back();
        if (TopologyInfo::dimension(type) != 3)
          continue;
        verts.clear();
        mesh->elements_get_attached_vertices( &e, 1, verts, junk, err ); MSQ_ERRRTN(err);
          
        for (unsigned s = 0; s < TopologyInfo::faces(type); ++s) {
          unsigned n;
          const unsigned *si = TopologyInfo::face_vertices( type, s, n );
          unsigned ns = 0;
          for (unsigned k = 0; k < n; ++k) {
            if (std::binary_search(surfaces[i]->vertices.begin(),
                                       surfaces[i]->vertices.end(),
                                       verts[si[k]]))
              ++ns;
          }
          if (ns >= 3) 
            sides.insert( std::pair<Mesh::ElementHandle,int>(e,s) );
        }
      }
    }
    
    std::vector< std::pair<Mesh::ElementHandle,int> > sstack;
    sstack.push_back( *sides.begin() );
    sides.erase(sides.begin());
    while (!sstack.empty()) {
      Mesh::ElementHandle e = sstack.back().first;
      int s = sstack.back().second;
      sstack.pop_back();
      
      verts.clear();
      mesh->elements_get_attached_vertices( &e, 1, verts, junk, err ); MSQ_ERRRTN(err);
      EntityTopology type;
      mesh->elements_get_topologies( &e, &type, 1, err ); MSQ_ERRRTN(err);
      unsigned n;
      const unsigned *si = TopologyInfo::face_vertices( type, s, n );
      
      // for each edge
      for (unsigned j = 0; j < n; ++j) {
        Mesh::VertexHandle v1 = verts[si[j]], v2 = verts[si[(j+1)%n]];
        vert_elems.clear();
        mesh->vertices_get_attached_elements( &v1, 1, vert_elems, junk, 
                                              err ); MSQ_ERRRTN(err);
        types.resize( vert_elems.size() );
        if (!vert_elems.empty()) {
          mesh->elements_get_topologies( arrptr(vert_elems), arrptr(types), 
                                 vert_elems.size(), err ); MSQ_ERRRTN(err);
        }
        while (!vert_elems.empty()) {
          Mesh::ElementHandle e2 = vert_elems.back();
          EntityTopology type2 = types.back();
          vert_elems.pop_back();
          types.pop_back();
          if (TopologyInfo::dimension(type) != 3)
            continue;
          verts2.clear();
          mesh->elements_get_attached_vertices( &e2, 1, verts2, junk, err ); MSQ_ERRRTN(err);
            // for each face
          for (unsigned s2 = 0; s2 < TopologyInfo::faces(type2); ++s2) {
            std::set< std::pair<Mesh::ElementHandle,int> >::iterator side;
            side = sides.find( std::pair<Mesh::ElementHandle,int>(e2,s2) );
            if (side == sides.end())
              continue;
          
            unsigned n2;
            const unsigned* si2 = TopologyInfo::face_vertices( type2, s2, n2 );
            unsigned idx;
            for (idx = 0; idx < n2; ++idx)
              if (verts2[si2[idx]] == v1)
                break;
            assert( idx < n2 );
            
            if (verts2[si2[(idx+1)%n2]] == v2 ||
                verts2[si2[(idx+n2-1)%n2]] == v2) {
              sstack.push_back( *side );
              sides.erase( side );
            }
          }
        }
      }
    }

    if (!sides.empty()) {
      MSQ_SETERR(err)("Surface mesh not a single, simply-connected patch", 
                      MsqError::INVALID_STATE);
      return;
    }
  }
}
    
  

static inline
bool operator<( const Mesh::EntityHandle h, 
                const DomainClassifier::DomainBlock& b )
  { return h < b.firstHandle; }

static inline
bool operator<( const DomainClassifier::DomainBlock& b,
                const Mesh::EntityHandle h )
  { return b.lastHandle < h; }

static inline
bool operator<( const DomainClassifier::DomainBlock& b,
			    const DomainClassifier::DomainBlock& c )
{ return b.lastHandle < c.firstHandle; }

MeshDomain* DomainClassifier::find_domain(
               Mesh::EntityHandle handle,
               const std::vector<DomainClassifier::DomainBlock>& list )
{
  std::vector<DomainClassifier::DomainBlock>::const_iterator i;
  i = std::lower_bound( list.begin(), list.end(), handle );
  return (i != list.end() && i->firstHandle <= handle) ? i->domain : NULL;
}

void DomainClassifier::snap_to( Mesh::EntityHandle entity_handle,
                                Vector3D &coordinate ) const
{
  if (const MeshDomain* dom = find_vertex_domain( entity_handle )) 
    dom->snap_to( entity_handle, coordinate );
}

void DomainClassifier::vertex_normal_at(Mesh::VertexHandle entity_handle,
                                  Vector3D &coordinate) const
{
  if (const MeshDomain* dom = find_vertex_domain( entity_handle )) 
    dom->vertex_normal_at( entity_handle, coordinate );
}

void DomainClassifier::element_normal_at(Mesh::ElementHandle entity_handle,
                                  Vector3D &coordinate) const
{
  if (const MeshDomain* dom = find_element_domain( entity_handle )) 
    dom->element_normal_at( entity_handle, coordinate );
}
                                  
void DomainClassifier::vertex_normal_at( const Mesh::VertexHandle* handles,
                                         Vector3D coordinates[],
                                         unsigned count,
                                         MsqError& err ) const
{
  for (unsigned i = 0; i < count; ++i) {
    const MeshDomain* dom = find_vertex_domain( handles[i] );
    if (!dom) {
      MSQ_SETERR(err)(MsqError::INVALID_ARG);
      return;
    }
    dom->vertex_normal_at( handles+i, coordinates+i, 1, err );
    MSQ_ERRRTN(err);
  }
}

void DomainClassifier::closest_point( Mesh::VertexHandle handle,
                                      const Vector3D& position,
                                      Vector3D& closest,
                                      Vector3D& normal,
                                      MsqError& err ) const
{
  if (const MeshDomain* dom = find_vertex_domain( handle )) 
    dom->closest_point( handle, position, closest, normal, err );
}

void DomainClassifier::domain_DoF( const Mesh::VertexHandle* handles,
                                   unsigned short* dof_array,
                                   size_t num_handles,
                                   MsqError& err ) const
{
  for (size_t i = 0; i < num_handles; ++i) {
    const MeshDomain* dom = find_vertex_domain( handles[i] );
    if (!dom) 
      dof_array[i] = 3;
    else {
      dom->domain_DoF( handles + i, dof_array + i, 1, err );
      MSQ_ERRRTN(err);
    }
  }
}

DomainClassifier::~DomainClassifier()
{
  if (deleteSubDomains)
    delete_all_sub_domains();
}

void DomainClassifier::delete_all_sub_domains()
{
    // get unique list of domains
  std::set<MeshDomain*> domains;
  std::vector<DomainBlock>::iterator i;
  for (i = vertexList.begin(); i != vertexList.end(); ++i)
    domains.insert( i->domain );
  for (i = elementList.begin(); i != elementList.end(); ++i)
    domains.insert( i->domain );
  std::set<MeshDomain*>::iterator j;
  for (j = domains.begin(); j != domains.end(); ++j)
    delete *j;
  clear();
}


} // namespace Mesquite
