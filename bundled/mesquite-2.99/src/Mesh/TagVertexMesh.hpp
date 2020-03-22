/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TagVertexMesh.hpp
 *  \brief Definition of Mesquite::TagVertexMesh class
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TAG_VERTEX_MESH_HPP
#define MSQ_TAG_VERTEX_MESH_HPP

#include "Mesquite.hpp"
#include "MeshDecorator.hpp"
#include "Instruction.hpp"

namespace MESQUITE_NS {

/**\brief Store alternate vertex coordinates in tags. 
 *
 * This class implements a decorator pattern for the Mesquite::Mesh
 * interface where alternate vertex coordinates are stored in a tag
 * on the original mesh.  The vertex coordinates are the same as that
 * of the decorated Mesh interface until they are set.  Once the coorindates
 * of a vertex are set using an instance of this class, the modified
 * vertex coordinates will be stored in a tag and returned from subsequent
 * queries of the vertex coordinates through this interface.
 *
 * The tag used to store alternate vertex coordinates is created and set
 * for all vertices when the coordinates of the any vertex are changed.
 * The tag type is a vector of three doubles.  
 *
 * Inserting an instance of this class into an InstructionQueue will result
 * in true vertex coordinates being copied into the alternate coordinate
 * values maintained by this class at that point in the instruction queue. 
 */
class MESQUITE_EXPORT TagVertexMesh : public MeshDecorator, public Instruction
{
  private:
  
    std::string tagName; //< Name of tag storing vertex coordinates
    TagHandle tagHandle;     //< Handle of tag storing vertex coordinates
    bool haveTagHandle;      //< True if tagHandle is set
    bool cleanUpTag;         //< If true, destroy tag in destructor
  
      /**\brief common code for constructor, set_mesh, and set_tag_name */
    void initialize( Mesh* mesh, std::string name, MsqError& );
      /**\brief copy real coordinate values into tag data */
    void copy_all_coordinates( MsqError& err );
      /**\brief if cleanUpTag, delete tag and clear handle */
    void check_remove_tag( MsqError& err );
  
  public:
  
    /** 
     *\param real_mesh  The mesh from which to aquire topology information
     *                  and vertex coordinates, and upon which to store
     *                  tags.
     *\param clean_up_tag_data If true, tag storing alternate vertex
     *                  coordinates will be removed when this object
     *                  is destroyed.
     *\param tag_name Name of tag in which to store alternate vertex coordinates.
     */
    TagVertexMesh( MsqError& err,
                   Mesh* real_mesh,
                   bool clean_up_tag_data = true,
                   std::string tag_name = "" );
    
      /** Destroy tag data for alternate coordinates if
       *  clean_up_tag_data is true. 
       */
    virtual ~TagVertexMesh();
    
    /**\brief Change the Mesh instance used as the real mesh. 
     *
     * Change the Mesh instance orignially specified in the 
     * constructor.
     * Note: Calling this function changes the handle space for
     *       mesh entities, invalidating any previous handle values, 
     *       iterators, etc. returned by the class instance.
     * Note: If clean_up_tag_data is true, calling this function
     *       will remove any stored alternate vertex coordinates 
     *        from the previous mesh.
     */
    void set_mesh( Mesh* real_mesh, MsqError& err );
    
    /**\brief Set tag cleanup behavior
     *
     * If true, class will remove any tag data storing alternate
     * vertex coordinates from the real mesh when a) the real Mesh
     * instance is changed or b) this object instance is destroted.
     */
    void should_clean_up_tag_data( bool value ) { cleanUpTag = value; }
    
    /**\brief Will tag storing alternate coordinates be destroyed. */
    bool will_clean_up_tag_data() const { return cleanUpTag; }
    
    /**\brief Get name of tag used to store alternate vertex coordinates. */
    std::string get_tag_name() const { return tagName; }

    /**\brief Set tag name used to store alternate vertex coordinates
     *
     * Change the tag name used to store alternate vertex coordinates.
     * Note:  Changing the tag name will result in the loss of any
     *        alternate vertex coordinates saved using the previous
     *        tag name.
     * Note:  If clean_up_tag_data is true, calling this function
     *        will result in the removal of the previous tag and
     *        any coordinates stored using that tag.
     *\param init  If the new tag already exists, any 
     *        coordinates stored in that tag will be used if this
     *        argument is false.  If this argument is true, the
     *        alternate coordinate values will be initialized to
     *        the true coordinate values in the real Mesh.
     */
    void set_tag_name( std::string name, MsqError& err );
  
    /**\brief clear all alternate vertex coordinate values
     *
     * Clear all alternate vertex coordinate values and 
     * revert to coordinates as stored in real mesh.
     */
    void clear( MsqError& err );
    
    

    virtual void vertices_get_coordinates( const VertexHandle vert_array[],
                                           MsqVertex* coordinates,
                                           size_t num_vtx,
                                           MsqError &err );

    virtual void vertex_set_coordinates( VertexHandle vertex,
                                         const Vector3D &coordinates,
                                         MsqError &err );
    

//***************  Tags  ***********

    virtual TagHandle tag_create( const std::string& tag_name,
                                  TagType type, unsigned length,
                                  const void* default_value,
                                  MsqError &err);

    virtual TagHandle tag_get( const std::string& name, 
                               MsqError& err );
    
//**************** Memory Management ****************

    virtual void release();

    
//**************** Instruction ****************

    virtual double loop_over_mesh( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError& err );

    virtual std::string get_name() const;
  
     //!\brief Called at start of instruction queue processing
     //!
     //! Do any preliminary global initialization, consistency checking,
     //! etc.  Default implementation does nothing.
    virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                   const Settings* settings,
                                   MsqError& err );
};


} // namespace Mesquite

#endif
