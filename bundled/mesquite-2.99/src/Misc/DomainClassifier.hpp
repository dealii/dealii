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

#ifndef MSQ_DOMAIN_CLASSIFIER_HPP
#define MSQ_DOMAIN_CLASSIFIER_HPP

/** \file DomainClassifier.hpp
 *  \brief 
 *  \author Jason Kraftcheck
 */

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

namespace MESQUITE_NS {

/**\brief Assign subsets of a mesh do different domains.
 *
 * Provide classification of mesh entities to domain(s).  For 
 * example, given a domain defined by multiple surfaces and
 * curves of intersection, associate the appropriate elements
 * and vertices to the appropriate geometric entities.
 */
class DomainClassifier : public MeshDomain
{
  public:
  
    /**\brief Check that classification maps to B-Rep topology */
    MESQUITE_EXPORT void test_valid_classification( Mesh* mesh, MsqError& err );
  
    /**\brief Classify mesh entities using tag values
     *
     * Given a list of MeshDomain instances, a unique
     * integer ID for each MeshDomain, and a tag name:
     * classify mesh entities by matching the tag value
     * on each mesh entity to the list of MeshDomain IDs.
     *
     *\param result   DomainClassifier to popupate
     *\param mesh     The Mesquite::Mesh instance
     *\param tag_name Tag containing integer domain ID for each mesh entity.
     *\param domain_array Array of MeshDomain instances.
     *\param id_array Array of integer MeshDomain IDs.
     *\param array_length Length of 'domain_array' and 'id_array'
     */
    static MESQUITE_EXPORT
    void classify_by_tag( DomainClassifier& result,
                          Mesh* mesh,
                          const char* tag_name,
                          MeshDomain** domain_array,
                          const int* id_array,
                          unsigned array_length,
                          MsqError& err );
  
    /**\brief Skin mesh and classify skin entities geometrically.
     *
     * Calculate the boundary of the mesh, and classify boundary
     * entities by calculating the distance from the mesh domain
     * to each vertex.
     *
     * NOTE:  Mesquite's limited MeshDomain callback interface
     * is not sufficient for robust geometric classification.  If
     * the mesh is not sufficiently refined, this method may produce
     * invalid results. 
     *
     * NOTE:  This method should not be used for surface meshes.
     * Everything in a surface mesh should be classified to a domain.
     * not just the skin.  Use classify_geometrically instead.
     *
     *\param result   DomainClassifier to popupate
     *\param mesh     The Mesquite::Mesh instance
     *\param tolerance Maximum distance a vertex may deviate from its domain.
     *\param domain_array Array of MeshDomain instances.
     *\param dimension_array Topological dimensiono of each MeshDomain
     *\param array_length Length of 'domain_array' and 'dimension_array'
     */
    MESQUITE_EXPORT static
    void classify_skin_geometrically( DomainClassifier& result,
                                      Mesh* mesh,
                                      double tolerance,
                                      MeshDomain** domain_array,
                                      const int* dimension_array,
                                      unsigned array_length,
                                      MsqError& err );

    /**\brief Classify all mesh entities geometrically.
     *
     * Classify entities by distance from vertex coordiantes to 
     * mesh domain instances.
     *
     * NOTE:  Mesquite's limited MeshDomain callback interface
     * is not sufficient for robust geometric classification.  If
     * the mesh is not sufficiently refined, this method may produce
     * invalid results. 
     *
     * NOTE:  This method is likely to fail for many domains for
     *        volume meshes.  Use classify_skin_geometrically for
     *        volume meshes.
     *
     *\param result   DomainClassifier to popupate
     *\param mesh     The Mesquite::Mesh instance
     *\param tolerance Maximum distance a vertex may deviate from its domain.
     *\param domain_array Array of MeshDomain instances.
     *\param dimension_array Topological dimensiono of each MeshDomain
     *\param array_length Length of 'domain_array' and 'dimension_array'
     */
   MESQUITE_EXPORT static
   void classify_geometrically( DomainClassifier& result,
                                Mesh* mesh,
                                double tolerance,
                                MeshDomain** domain_array,
                                const int* dimension_array,
                                unsigned array_length,
                                MsqError& err );
  
    struct DomainSet {
      MESQUITE_EXPORT DomainSet( MeshDomain* dom ) : domain(dom) {}
      MESQUITE_EXPORT DomainSet() : domain(0) {}
      MeshDomain* domain;
	  MESQUITE_EXPORT void set_vertices( const std::vector<Mesh::VertexHandle>& verts ) 
	    { vertices = verts; }
	  MESQUITE_EXPORT void set_elements( const std::vector<Mesh::ElementHandle>& elems ) 
	    { elements = elems; }
	  MESQUITE_EXPORT void get_vertices( std::vector<Mesh::VertexHandle>& verts ) const
	    { verts = vertices; }
	  MESQUITE_EXPORT void get_elements( std::vector<Mesh::ElementHandle>& elems ) const
	    { elems = elements; }
      std::vector<Mesh::VertexHandle> vertices;
      std::vector<Mesh::ElementHandle> elements;
    };
  
    /**\brief Specify classification explicitly for each entity.
     *
     * For each mesh element and vertex that is classified to
     * a domain, specify that classification explicitly.
     *
     *\param result           DomainClassifier to popupate
     *\param mesh             The Mesquite::Mesh instance
     *\param domain_set_array Array of DomainClassifier::DomainSet structs
     *                        specifying for each MeshDomain: the
     *                        domain instance and the elements and
     *                        vertices associated with the domain.
     *\param array_length     Length of 'domain_set_array'
     */
    MESQUITE_EXPORT static
    void classify_by_handle( DomainClassifier& result,
                             Mesh* mesh,
                             DomainSet* domain_set_array,
                             unsigned array_length,
                             MsqError& err );
 
    MESQUITE_EXPORT DomainClassifier() : deleteSubDomains(false)
      {}
    
    MESQUITE_EXPORT virtual ~DomainClassifier();
    
    MESQUITE_EXPORT
    virtual void snap_to( Mesh::VertexHandle entity_handle,
                          Vector3D &coordinate) const;
    
    MESQUITE_EXPORT
    virtual void vertex_normal_at( Mesh::VertexHandle entity_handle,
                                   Vector3D &coordinate) const;
    MESQUITE_EXPORT
    virtual void element_normal_at( Mesh::ElementHandle entity_handle,
                                    Vector3D &coordinate) const;
                          
    MESQUITE_EXPORT
    virtual void vertex_normal_at( const Mesh::VertexHandle* handles,
                                   Vector3D coordinates[],
                                   unsigned count,
                                   MsqError& err ) const;
                            
    MESQUITE_EXPORT
    virtual void closest_point( Mesh::VertexHandle handle,
                                const Vector3D& position,
                                Vector3D& closest,
                                Vector3D& normal,
                                MsqError& err ) const;
                                
    MESQUITE_EXPORT
    virtual void domain_DoF( const Mesh::VertexHandle* handle_array,
                             unsigned short* dof_array,
                             size_t num_handles,
                             MsqError& err ) const;
    
      /**\brief Clear all data, including MeshDomain list */
    MESQUITE_EXPORT
    void clear() {
      vertexList.clear();
      elementList.clear();
    }
  
    struct DomainBlock {
      Mesh::EntityHandle firstHandle;
      Mesh::EntityHandle lastHandle;
      MeshDomain* domain;
    };
                                      
 	MESQUITE_EXPORT
    MeshDomain* find_vertex_domain( Mesh::VertexHandle vertex ) const
      { return find_domain( vertex, vertexList ); }
	MESQUITE_EXPORT
    MeshDomain* find_element_domain( Mesh::ElementHandle element ) const
      { return find_domain( element, elementList ); }
      
    MESQUITE_EXPORT
    void delete_sub_domains(bool yesno)
      { deleteSubDomains = yesno; }

    MESQUITE_EXPORT
    void delete_all_sub_domains();
    
  private:
    bool deleteSubDomains;
    
    static MeshDomain* find_domain( Mesh::EntityHandle handle, 
                                 const std::vector<DomainBlock>& list );
    
    std::vector<DomainBlock> vertexList;
    std::vector<DomainBlock> elementList;
};

} // namespace Mesquite

#endif // MSQ_DOMAIN_CLASSIFIER_HPP
