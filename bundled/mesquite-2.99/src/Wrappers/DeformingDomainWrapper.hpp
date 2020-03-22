/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file DeformingDomainWrapper.hpp
 *  \brief Define DeformingDomainWrapper class
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_DEFORMING_DOMAIN_WRAPPER_HPP
#define MSQ_DEFORMING_DOMAIN_WRAPPER_HPP

#include "Wrapper.hpp"
#include "MeshInterface.hpp"
#include <string>

namespace MESQUITE_NS {

class CurveDomain;

/**\brief Smooth mesh with deforming domain (deforming geometry)
 *
 * This wrapper stores the initial mesh on/in an undeformed domain
 * and uses the characteristics of that mesh to guide the optimization
 * of the mesh after deformation of the domain.  
 *
 * The caller must ensure that the wrapper has the initial mesh 
 * description.  The simplest way to do this is to call \c store_initial_mesh
 * prior to the deformation of the geometric domain. 
 *
 * Application using this wrapper should preform the following steps.
 * -# Change any wrapper settings
 * -# Call \c store_initial_mesh with a mesh object containing all vertices in the mesh
 * -# If using \c DeformingCurveSmoother, call DeformingCurveSmoother::store_initial_mesh
 *      separately for each geometric curve in the mesh.
 * -# Optionally export mesh to file if \c Mesh implementation can save tag data
 * -# Deform domain
 * -# Optionally restore mesh from file
 * -# Move any vertices bound to geometric points/vertices to their new location
 * -# If using \c DeformingCurveSmoother, call \c smooth_curve separately for each geometric curve
 * -# If doing free smooth call \c run_instructions for entire mesh.  Otherwise
 *       call run_instructions separately for each geometric surface and each
 *       volume.
 *
 *\NOTE This algorithm uses non-barrier metrics which makes it possible, if
 *      unlikely, that the result mesh will contain more/different inverted
 *      elements than the input mesh.
 *
 *\NOTE Mesquite does not do edge/curve smoothing.  The caller may either
 *      set up a "free" smooth where vertices on geometric curves are 
 *      smoothed as a part of teh surface or volume optimization or simply
 *      redistribute the vertices along the curves outside of mesquite.
 *      The DeformingCurveSmoother class is provided below to facilitate
 *      redistributing nodes on a curve such that fractional arc length
 *      is preserved for edges.  However, this it is as a stand-alone utility 
 *      that must be invoked independently for each curve.
 */
class DeformingDomainWrapper : public Wrapper
{
public:

  MESQUITE_EXPORT
  DeformingDomainWrapper();
  
  MESQUITE_EXPORT virtual
  ~DeformingDomainWrapper();
  
  /**\brief Specify tag used to get/set cached initial mesh vertex coordinates
   *
   * This optimizer uses the initial mesh (prior to deformation of the 
   * domain) as a reference when smoothing the mesh on the deformed domain.
   * The initial mesh is stored/accessed as alternate vertex coordinates
   * stored in a Tag on the mesh.  This function can be used to specify
   * the tag used to store those initial coordinates.  If the named tag
   * already exists then it must have a type of DOUBLE and a size of three values.  
   *
   * If the application already has the initial coordinates stored in some
   * tag, then calling this function with the tag handle is sufficient to
   * to set up the optmizer.  Otherwise \c store_initial_mesh can be called
   * to copy the initial vertex coordinates into a tag.  If both this function
   * is called to set the tag handle and \c store_initial_mesh is called, then
   * \c store_initial_mesh will copy the vertices into the specified tag.
   * Otherwise \c store_initial_mesh, if called, will create a new tag.
   */
  inline
  void set_initial_coord_tag( const std::string& tag_name )
    { initVertexCoords = tag_name; }
  
  /**\brief Get tag used to store vertex coordinates of initial mesh
   *
   * Return the tag passed to \c set_initial_coord_tag or if
   * \c store_initial_mesh has been called, the the tag created by
   * \c store_initial_mesh .
   */
  inline
  std::string get_initial_coord_tag() const
    { return initVertexCoords; }
  
  /**\brief Store the initial mesh before deformation as a reference
   *
   * Store the initial mesh for use as a reference when smoothing
   * mesh on deformed geometry.
   *\param mesh The mesh instance to be operated on.  Initial mesh
   *      data is stored by copying vertex coordinates to 
   *\Note This must be called prior to the optimization with the
   *      mesh in its undeformed state.  
   *\Note This function need not be called if the application cal
   *      provide access to the coordinates of vertices of the initial
   *      mesh through the tag mechanism already.  
   *      See \c set_initial_coord_tag
   */
  MESQUITE_EXPORT
  void store_initial_mesh( Mesh* mesh, MsqError& err );
  
  /**\brief Mesh characteristics to attempt to preserve */
  enum MeshCharacteristic { 
    SHAPE,              //!< Attempt to preserve element shape (TShapeNB1 metric)
    SHAPE_SIZE,         //!< Attempt to preserve element shape and size (TShapeSize2DNB1 & TShapeSize3DNB1)
    SHAPE_SIZE_ORIENT   //!< Attempt to preserve element shape, size, and orientation (TShapeSizeOrientNB1)
  }; 
  
  /**\brief Specify which mesh characteristics to preserve 
   *
   * It typically works the best to preserve whichever quantities can
   * be preserved.  For example, if the overall area/volume of the domain
   * didn't change much, attempt to preserve size.  If the deformation of
   * the domain doesn't require rotation of the elements relative to the
   * global coordinate system, then preserve orientation.
   *
   * The default is \c SHAPE.
   */
  inline
  void set_mesh_characteristic( MeshCharacteristic t )
    { metricType = t; }
  
  inline
  MeshCharacteristic get_mesh_characteristic() const
    { return metricType; }
  
  /**\brief Check if vertex culling will be used */
  inline bool is_culling_enabled() const
    { return doCulling; }
  
  /**\brief Enable vertex culling.
   * 
   * Culling will be enabled by default.
   */
  inline void enable_culling( bool yesno )
    { doCulling = yesno; }

  /**\brief Specify timeout after which optimization will exit */
  inline
  void set_cpu_time_limit( double seconds )
    { cpuTime = seconds; }

  /**\brief No timeout */
  inline
  void clear_cpu_time_limit()
    { cpuTime = 0.0; }
  
  /**\brief Get timeout after which optimization will exit
   *\return false if no timeout specified, true otherwise */
  inline 
  bool get_cpu_time_limit( double& seconds ) const 
    { seconds = cpuTime; return cpuTime > 0.0; }

  /**\brief Specify factor by which to minimum distance a vertex must 
   *        move in an iteration to avoid termination of the optimization.
   *
   * Value must be greater than zero, and should be less 1.0
   * The optimizer will calculated a representative edge length
   * from the mesh.  It will terminate the optimization when an
   * iteration moves no vertex more than this factor times the
   * characteristic edge length.
   */
  MESQUITE_EXPORT
  void set_vertex_movement_limit_factor( double f );
  
  inline 
  double get_vertex_movement_limit_factor( ) const 
    { return movementFactor; }  

protected:

  MESQUITE_EXPORT
  void run_wrapper( MeshDomainAssoc* mesh_and_domain,
                    ParallelMesh* pmesh,
                    Settings* settings,
                    QualityAssessor* qa,
                    MsqError& err );

  MESQUITE_EXPORT
  void move_to_domain( Mesh* mesh, MeshDomain* geom, MsqError& err );
                    

private:

  MeshCharacteristic metricType;
  std::string initVertexCoords; //!< name of tag used to store coordinates of initial mesh
  bool doCulling;
  double cpuTime, movementFactor;
};

/**\brief Utility to do curve smoothing
 *
 * This class implements a simple utility to redistribute vertices
 * on a curve (1D domain).  It can operate in two modes.  If the 
 * charactristic is \c EQUAL, it will redistribute the vertices such
 * that the arc length between each vertex is equal. In the 
 * \c PROPORTIONAL mode the application must call \c store_initial_mesh
 * for each curve before the domain is deformed.  \c store_initial_mesh
 * will record the fractional arc length between vertices so that
 * \c smooth_curve can redistribute the vertices such that the space
 * between the vertices is the same fractional arc length.
 */
class DeformingCurveSmoother 
{
public:

  MESQUITE_EXPORT
  DeformingCurveSmoother();
  
  MESQUITE_EXPORT
  ~DeformingCurveSmoother();
  
  /**\brief Mesh characteristics to attempt to preserve */
  enum Scheme { 
     EQUAL,       //!< space curve vertices equally
     PROPORTIONAL //!< preserve relative spacing from initial mesh
  };
  
  /**\brief Specify tag used to get/set cached initial mesh data
   *
   * This optimizer uses the initial mesh (prior to deformation of the 
   * domain) as a reference when smoothing the mesh on the deformed domain.
   * This function can be used to store initial mesh data.  If the tag 
   * already exists then t must have a type of DOUBLE and a size of one value.
   * the tag used to store those initial coordinates.  It must have a type
   * of DOUBLE and a size of one value.  
   */
  inline
  void set_initial_fraction_tag( const std::string& tag_name )
    { initFractTag = tag_name; }
  
  /**\brief Get tag used to store initial mesh characteristics
   *
   * Return the tag handle passed to \c set_initial_fraction_tag or if
   * \c store_initial_mesh, the the tag created by in.
   */
  inline
  std::string get_initial_fraction_tag() const
    { return initFractTag; }
  
  /**\brief Specify which mesh characteristics to preserve */
  inline
  void set_mesh_characteristic( Scheme t )
    { metricType = t; }
  
  inline
  Scheme get_mesh_characteristic() const
    { return metricType; }
  
  /**\brief Record relative edge length from initial mesh
   *
   * Record the necessary characteristics of the initial mesh required
   * for a later all to \c smooth_curve using the \c PROPORTIONAL mode.
   *\param mesh_instance  This mesh instance need not contain only the
   *               vertices of this curve, but it must preserve tag data
   *               on the curve vertices for use in \c smooth_curve.
   *\param vertex_array   The array of handles corresponding to vertices
   *               on the curve and its end points.  Vertices must be 
   *               ordered and in the forward direction along the curve,
   *               where the forward direction is whichever direction
   *               the \c CurveDomain instance considers a positive arc
   *               length for the \c position_from_length call.  The array
   *               must also contain both the start and end vertex for 
   *               the curve, and if the curve is closed, the single start/end
   *               vertex must occur at both the start and end of the list.
   *\param vertex_array_length The number of handles in \c vertex_array.
   *\param geometry The curve or 1D domain geometry.
   */
  MESQUITE_EXPORT
  void store_initial_mesh( Mesh* mesh_instance,
                           const Mesh::VertexHandle* vertex_array,
                           int vertex_array_length,
                           CurveDomain* geometry,
                           MsqError& err );
  
  /**\brief Redistribute vertices along a curve or 1D domain.
   *
   * Smooth the curve mesh by redistributing the vertices using one of
   * the schemes defined in \c Scheme.
   *
   *\param mesh_instance  This mesh instance need not contain only the
   *               vertices of this curve, but if using the 
   *               \c PROPORTIONAL scheme it must have the tag data
   *               defined by the earlier call to \c store_initial_mesh.
   *\param vertex_array   The array of handles corresponding to vertices
   *               on the curve and its end points.  Vertices must be 
   *               ordered and in the forward direction along the curve,
   *               where the forward direction is whichever direction
   *               the \c CurveDomain instance considers a positive arc
   *               length for the \c position_from_length call.  The array
   *               must also contain both the start and end vertex for 
   *               the curve, and if the curve is closed, the single start/end
   *               vertex must occur at both the start and end of the list.
   *               The first and last vertex in the array will *not* be moved.
   *               If the geometric points (0D domains) they are constrained
   *               to lie on have moved, the caller must move them to those
   *               new locations before calling this function.
   *\param vertex_array_length The number of handles in \c vertex_array.
   *\param geometry The curve or 1D domain geometry.
   *\param type    The scheme to use.  If using \c PROPORTIONAL, 
   *               then \c store_initial_mesh must be called for the same
   *               vertices and mesh instance before any deformation.
   */
  MESQUITE_EXPORT
  void smooth_curve( Mesh* mesh_instance,
                     const Mesh::VertexHandle* vertex_array,
                     int vertex_array_length,
                     CurveDomain* geometry,
                     Scheme type,
                     MsqError& err );

private:
  Scheme metricType;
  std::string initFractTag; //!< name of tag used to store arc length fractions
  TagHandle get_tag(Mesh* mesh, MsqError& err);
};


} // namespace MESQUITE_NS

#endif
