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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov,
    kraftche@cae.wisc.edu
   
  ***************************************************************** */
/*!
  \file   MeshInterface.hpp
  \brief  This file contains the Mesquite mesh interface.
          Many users will want to implement a concrete class derived from
          the MeshInterface class to access their mesh.


  \author Darryl Melander
  \author Thomas Leurent
  \date   2003-04-17
*/
#ifndef MESQUITE_INTERFACE_HPP
#define MESQUITE_INTERFACE_HPP

#include "Mesquite.hpp"
#include "TopologyInfo.hpp"

#include "MsqError.hpp"
#include "MsqVertex.hpp"

#include <vector>
#include <cstddef>
#include <string>
#include <iostream>

namespace MESQUITE_NS
{
  class EntityIterator;
  class MsqError;
  class MsqVertex;
  class Vector3D;
  class SphericalDomain;
  typedef EntityIterator VertexIterator;
  typedef EntityIterator ElementIterator;
    
    /** Type used to refer to a tag defintion */
  typedef void* TagHandle;

  inline size_t vertices_in_topology(EntityTopology);

  /*! \class Mesh
     \brief  A Mesquite::Mesh is a collection of mesh elements which are
     composed of mesh vertices.  Intermediate objects are not accessible
     through this interface (where intermediate objects include things
     like the faces of a hex, or an element's edges).
  */
  class MESQUITE_EXPORT Mesh
  {
  public:
//************ Type Definitions **************
      //! Opaque EntityHandle type and tag type.
    typedef void* EntityHandle;
    
    
      // We typedef specific types of EntityHandles just
      // to make it clear what kind of entity is to be
      // returned or used as a function parameter, but the
      // different handle types are not actually distinct.
    typedef EntityHandle VertexHandle;
    typedef EntityHandle ElementHandle;

//************ Operations on entire mesh ****************
      //! Returns whether this mesh lies in a 2D or 3D coordinate system.
    virtual int get_geometric_dimension(MsqError &err) = 0;
    
    /** \brief Get all elements in mesh
     *
     * Get the handles of every element in the active mesh.
     */ 
    virtual void get_all_elements( std::vector<ElementHandle>& elements,
                                   MsqError& err ) = 0;
    
    /** \brief Get all vertices in mesh
     *
     * Get the handles of every vertex in the active mesh
     */
    virtual void get_all_vertices( std::vector<VertexHandle>& vertices,
                                   MsqError& err ) = 0;
    
      //! Returns a pointer to an iterator that iterates over the
      //! set of all vertices in this mesh.  The calling code should
      //! delete the returned iterator when it is finished with it.
      //! If vertices are added or removed from the Mesh after obtaining
      //! an iterator, the behavior of that iterator is undefined.
//    virtual VertexIterator* vertex_iterator(MsqError &err) = 0;
    
      //! Returns a pointer to an iterator that iterates over the
      //! set of all top-level elements in this mesh.  The calling code should
      //! delete the returned iterator when it is finished with it.
      //! If elements are added or removed from the Mesh after obtaining
      //! an iterator, the behavior of that iterator is undefined.
//    virtual ElementIterator* element_iterator(MsqError &err) = 0;

//************ Vertex Properties ********************
      //! Returns true or false, indicating whether the vertex
      //! is allowed to be repositioned.  True indicates that the vertex
      //! is fixed and cannot be moved.  Note that this is a read-only
      //! property; this flag can't be modified by users of the
      //! Mesquite::Mesh interface.
    virtual void vertices_get_fixed_flag( const VertexHandle vert_array[], 
                                          std::vector<bool>& fixed_flag_array,
                                          size_t num_vtx, 
                                          MsqError &err ) = 0;

      //! Returns true or false, indicating whether the vertex
      //! is a higher-order node that should be slaved to the logical
      //! mid-point of the element side it lies on or not, respectively.  
      //!
      //! Note: This function will never be called unless this behavior is
      //! requested by calling:
      //! InstructionQueue::set_slaved_ho_node_mode( Settings::SLAVE_FLAG )
    virtual void vertices_get_slaved_flag( const VertexHandle vert_array[], 
                                           std::vector<bool>& slaved_flag_array,
                                           size_t num_vtx, 
                                           MsqError &err ) = 0;

      //! Get/set location of a vertex
    virtual void vertices_get_coordinates( const VertexHandle vert_array[],
                                           MsqVertex* coordinates,
                                           size_t num_vtx,
                                           MsqError &err ) = 0;
    virtual void vertex_set_coordinates( VertexHandle vertex,
                                         const Vector3D &coordinates,
                                         MsqError &err ) = 0;
    
      //! Each vertex has a byte-sized flag that can be used to store
      //! flags.  This byte's value is neither set nor used by the mesh
      //! implementation.  It is intended to be used by Mesquite algorithms.
      //! Until a vertex's byte has been explicitly set, its value is 0.
    virtual void vertex_set_byte( VertexHandle vertex,
                                  unsigned char byte, 
                                  MsqError &err) = 0;
    virtual void vertices_set_byte( const VertexHandle *vert_array,
                                    const unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err ) = 0;
    
      //! Retrieve the byte value for the specified vertex or vertices.
      //! The byte value is 0 if it has not yet been set via one of the
      //! *_set_byte() functions.
    virtual void vertex_get_byte( const VertexHandle vertex,
                                  unsigned char *byte, 
                                  MsqError &err ) = 0;
    virtual void vertices_get_byte( const VertexHandle *vertex,
                                    unsigned char *byte_array,
                                    size_t array_size, 
                                    MsqError &err ) = 0;
    
//**************** Vertex Topology *****************    
      /** \brief get elements adjacent to vertices
       *
       * Get adjacency data for vertices
       *
       *\param vertex_array    Array of vertex handles specifying the
       *                       list of vertices to retrieve adjacency
       *                       data for.
       *\param num_vertex      Number of vertex handles in #vertex_array
       *\param elements     The array in which to place the handles of
       *                       elements adjacent to the input vertices.
       *\param offsets    For each vertex in #vertex_array, the
       *                       value in the corresponding position in this
       *                       array is the index into #elem_array at
       *                       which the adjacency list begins for that
       *                       vertex.
       */
    virtual void vertices_get_attached_elements( 
                         const VertexHandle* vertex_array,
                         size_t num_vertex,
                         std::vector<ElementHandle>& elements,
                         std::vector<size_t>& offsets,
                         MsqError& err ) = 0;
    
//*************** Element Topology *************
    
      /** \brief Get element connectivity
       *
       * Get the connectivity (ordered list of vertex handles) for
       * each element in the input array.
       *
       *\param elem_handles  The array of element handles for which to
       *                     retrieve the connectivity list.
       *\param num_elems     The length of #elem_handles
       *\param vert_handles  Array in which to place the vertex handles
       *                     in each elements connectivity.
       *\param offsets       For each element in #elem_handles, the
       *                     value in the same position in this array
       *                     is the index into #vert_handles at which
       *                     the connectivity list for that element begins.
       */
    virtual void elements_get_attached_vertices(
                                   const ElementHandle *elem_handles,
                                   size_t num_elems,
                                   std::vector<VertexHandle>& vert_handles,
                                   std::vector<size_t>& offsets, 
                                   MsqError &err) = 0;
    
    
      //! Returns the topologies of the given entities.  The "entity_topologies"
      //! array must be at least "num_elements" in size.
    virtual void elements_get_topologies(const ElementHandle *element_handle_array,
                                         EntityTopology *element_topologies,
                                         size_t num_elements, MsqError &err) = 0;

    
//***************  Tags  ***********
    
      /** The type of a tag */
    enum TagType { BYTE, BOOL, INT, DOUBLE, HANDLE };

      /** \brief Create a tag
       *
       * Create a user-defined data type that can be attached
       * to any element or vertex in the mesh.  For an opaque or
       * undefined type, use type=BYTE and length=sizeof(..).
       *
       * \param tag_name  A unique name for the data object
       * \param type      The type of the data
       * \param length    Number of values per entity (1->scalar, >1 ->vector)
       * \param default_value Default value to assign to all entities - may be NULL
       * \return - Handle for tag definition 
       */
    virtual TagHandle tag_create( const std::string& tag_name,
                                  TagType type, unsigned length,
                                  const void* default_value,
                                  MsqError &err) = 0;
    
      /** \brief Remove a tag and all corresponding data
       *
       * Delete a tag.
       */
    virtual void tag_delete( TagHandle handle, MsqError& err ) = 0;
    
    
      /** \brief Get handle for existing tag, by name. 
        *
        * Check for the existance of a tag given it's name and
        * if it exists return a handle for it.  If the specified
        * tag does not exist, zero should be returned WITHOUT 
        * flagging an error.
        */
    virtual TagHandle tag_get( const std::string& name, 
                               MsqError& err ) = 0;
     
      /** \brief Get properites of tag
       *
       * Get data type and number of values per entity for tag.
       * \param handle     Tag to get properties of.
       * \param name_out   Passed back tag name.
       * \param type_out   Passed back tag type.
       * \param length_out Passed back number of values per entity.
       */
    virtual void tag_properties( TagHandle handle,
                                 std::string& name_out,
                                 TagType& type_out,
                                 unsigned& length_out,
                                 MsqError& err ) = 0;
    
      /** \brief Set tag values on elements
       * 
       * Set the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       const void* tag_data,
                                       MsqError& err ) = 0;

      /** \brief Set tag values on vertices
       * 
       * Set the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of node_array
       * \param node_array Array of vertices for which to set the tag value.
       * \param tag_data   Tag data for each element, contiguous in memory.
       *                   This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_set_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       const void* tag_data,
                                       MsqError& err ) = 0;
    
    
      /** \brief Get tag values on elements
       * 
       * Get the value of a tag for a list of mesh elements.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of elements for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_element_data( TagHandle handle,
                                       size_t num_elems,
                                       const ElementHandle* elem_array,
                                       void* tag_data,
                                       MsqError& err ) = 0;
    
      /** \brief Get tag values on vertices.
       * 
       * Get the value of a tag for a list of mesh vertices.
       * \param handle     The tag 
       * \param num_elems  Length of elem_array
       * \param elem_array Array of vertices for which to get the tag value.
       * \param tag_data   Return buffer in which to copy tag data, contiguous 
       *                   in memory.  This data is expected to be 
       *                   num_elems*tag_length*sizeof(tag_type) bytes.
       */
    virtual void tag_get_vertex_data ( TagHandle handle,
                                       size_t num_elems,
                                       const VertexHandle* node_array,
                                       void* tag_data,
                                       MsqError& err ) = 0;

    
    
//**************** Memory Management ****************
      //! Tells the mesh that the client is finished with a given
      //! entity handle.  
    virtual void release_entity_handles( const EntityHandle *handle_array,
                                         size_t num_handles, 
                                         MsqError &err) = 0;
    
      //! Instead of deleting a Mesh when you think you are done,
      //! call release().  In simple cases, the implementation could
      //! just call the destructor.  More sophisticated implementations
      //! may want to keep the Mesh object to live longer than Mesquite
      //! is using it.
    virtual void release() = 0;
    
    virtual ~Mesh()
      {}
  protected:
  };
  
  /*! \class EntityIterator 
  \brief   Iterates through a set of entities.  An EntityIterator 
  is typically obtained via Mesh::vertex_iterator() or
  Mesh::element_iterator().  An iterator obtained in this
  way iterates over the set of all vertices/elements in
  the Mesh from which the iterator was obtained. */
  class EntityIterator
  {
  public:
    virtual ~EntityIterator()
      {}

      //! Moves the iterator back to the first
      //! entity in the list.
    virtual void restart() = 0;
    
      //! *iterator.  Return the handle currently
      //! being pointed at by the iterator.
    virtual Mesh::EntityHandle operator*() const = 0;
    
      //! ++iterator
    virtual void operator++() = 0;
    
      //! Returns false until the iterator has
      //! been advanced PAST the last entity.
      //! Once is_at_end() returns true, *iterator
      //! returns 0.
    virtual bool is_at_end() const = 0;
  };
  

  /*! \class MeshDomain
      The MeshDomain class provides geometrical information concerning the Mesh.
      It is called during surface meshes optimization to figure out the surface normal,
      how to snap vertices back to the surface, etc... . 
    */
  class MESQUITE_EXPORT MeshDomain
  {
  public:
    virtual ~MeshDomain()
      {}
      
      //! Modifies "coordinate" so that it lies on the
      //! domain to which "entity_handle" is constrained.
      //! The handle determines the domain.  The coordinate
      //! is the proposed new position on that domain.
    virtual void snap_to(Mesh::VertexHandle entity_handle,
                         Vector3D &coordinate) const = 0;
    
      //! Returns the normal of the domain to which
      //! "entity_handle" is constrained.  For non-planar surfaces,
      //! the normal is calculated at the point on the domain that
      //! is closest to the passed in value of "coordinate".  If the
      //! domain does not have a normal, or the normal cannot
      //! be determined, "coordinate" is set to (0,0,0).  Otherwise,
      //! "coordinate" is set to the domain's normal at the
      //! appropriate point.
      //! In summary, the handle determines the domain.  The coordinate
      //! determines the point of interest on that domain.
      //!
      //! User should see also PatchData::get_domain_normal_at_vertex and
      //! PatchData::get_domain_normal_at_element .
    virtual void vertex_normal_at(Mesh::VertexHandle entity_handle,
                                  Vector3D &coordinate) const = 0;
    virtual void element_normal_at(Mesh::ElementHandle entity_handle,
                                  Vector3D &coordinate) const = 0;
                          
      /**\brief evaluate surface normals
       *
       * Returns normals for a domain.
       *
       *\param handles       The domain evaluated is the one in which
       *                     this mesh entity is constrained.
       *\param coordinates   As input, a list of positions at which to
       *                     evaluate the domain.  As output, the resulting
       *                     domain normals.
       *\param count         The length of the coordinates array.
       */
    virtual void vertex_normal_at( const Mesh::VertexHandle* handles,
                                   Vector3D coordinates[],
                                   unsigned count,
                                   MsqError& err ) const = 0;
                            
      /**\brief evaluate closest point and normal
       *
       * Given a position in space, return the closest 
       * position in the domain and the domain normal
       * at that point.
       *
       *\param entity_handle Evaluate the subset of the domain contianing
       *                     this entity
       *\param position      Input position for which to evaluate
       *\param closest       Closest position in the domain.
       *\param normal        Domain normal at the location of 'closest'
       */
    virtual void closest_point( Mesh::VertexHandle handle,
                                const Vector3D& position,
                                Vector3D& closest,
                                Vector3D& normal,
                                MsqError& err ) const = 0;
                                
      /**\brief Get degrees of freedom in vertex movement.
       *
       * Given a vertex, return how the domain constrains the
       * location of that vertex as the number of degrees of
       * freedom in the motion of the vertex.  If the domain
       * is a geometric domain, the degrees of freedom for a
       * vertex is the dimension of the geometric entity the
       * vertex is constrained to lie on (e.g. point = 0, curve = 1,
       * surface = 2, volume = 3.)
       */
    virtual void domain_DoF( const Mesh::EntityHandle* handle_array,
                             unsigned short* dof_array,
                             size_t num_handles,
                             MsqError& err ) const = 0;
                             
       
  };

  /*! \class MeshDomainAssoc
      The MeshDomainAssoc class provides an association of a Mesh instance
      with a MeshDomain instance.  The mesh is checked to verify that
      it is compatibile with the associated MeshDomain.  If the two are
      not compatible, the MeshDomainAssoc instace is not created. 
    */
  class MESQUITE_EXPORT MeshDomainAssoc
  {
  public:

      /**\brief Constructor
       *\param mesh                       The mesh instance being associated
       * param domain                     The domain being associated
       * param full_compatibility_check   Controls how many vertices will be checked for 
       *                                  compatibility with the associated domain.
       *                                  When false, only the first vertex of the mesh
       *                                  is checked for compatibility.  When true, all
       *                                  vertices of the mesh are checked.
       * param proceed                    Controls what Mesquite will do if the compatibility
       *                                  check fails.  When false, mesquite terminates i
       *                                  execution.  When true, execution continues.
       * param skip_compatibility_checki  when true, does not perform the compatibility check. 
       *                                  When false, the check is performed.
       */
    MeshDomainAssoc(Mesquite::Mesh* mesh, 
                    Mesquite::MeshDomain* domain, 
                    bool full_compatibility_check=false,
                    bool proceed=false,
                    bool skip_compatibility_check=false)
       : mMesh(mesh), mMeshDomain(domain), mesh_and_domain_are_compatible(false)
    {
        // check for real instance.  If either value is NULL then it's just an 
        // instance created to facilitate passing of just a mesh or domain
        // also, check if skipping the compatibility check was requested
      if (mesh && domain && !skip_compatibility_check)
      {
        MsqError err;
        double tolerance = 1.0e-3;

        std::vector<Mesh::VertexHandle> vert_handles;
        mMesh->get_all_vertices( vert_handles, err );     

        MsqVertex mesh_vertex, domain_vertex;
        Vector3D normal;

        double distance; 
        std::vector<int>::size_type i, times_to_loop;
        if (full_compatibility_check)
          times_to_loop = vert_handles.size();
        else
          times_to_loop = 1;
        mesh_and_domain_are_compatible = true;
        for (i = 0; i < times_to_loop; ++i) 
        {
          mMesh->vertices_get_coordinates(&vert_handles[i], 
                                          &mesh_vertex,
                                          1,
                                          err);     
          mMeshDomain->closest_point( vert_handles[i],
                                      Vector3D(mesh_vertex),
                                      domain_vertex,
                                      normal,
                                      err ); 

        distance = Vector3D::distance_between(mesh_vertex, domain_vertex);
        if ( distance > tolerance )
        {
          mesh_and_domain_are_compatible = false;
          std::cout << "Warning: Mesh and Domain are not compatibile" << std::endl;
          if (!proceed)
          {
            std::cout << "Terminating due to Mesh/Domain incompatibility" << std::endl;
            throw "Terminating due to Mesh/Domain incompatibility";
          }
          break;   // exits for loop when not compatbile but should not terminate
        }
      }
    }
  }


    ~MeshDomainAssoc() {};
   
  
       /**\brief get associated mesh
       *
       * Return the mesh associated with this instance.
       */
    Mesquite::Mesh* get_mesh()
    {
      return mMesh;
    };
   

      /**\brief get associated domain
       *
       * Return the domain associated with this instance.
       */
    Mesquite::MeshDomain* get_domain()
    {
      return mMeshDomain;
    };


  private:
    Mesquite::Mesh* mMesh;
    Mesquite::MeshDomain* mMeshDomain;
    bool mesh_and_domain_are_compatible;

  public:

    bool are_compatible()
    {
      return mesh_and_domain_are_compatible;
    };

  };
}

inline size_t Mesquite::vertices_in_topology(Mesquite::EntityTopology topo)
{
  return TopologyInfo::corners( topo );
}

#endif
