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

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file Settings.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_SETTINGS_HPP
#define MSQ_SETTINGS_HPP

#include "Mesquite.hpp"
#include <stdlib.h> // for size_t

namespace MESQUITE_NS {

class MappingFunction;
class MappingFunction2D;
class MappingFunction3D;
struct SettingData;

class MESQUITE_EXPORT Settings {
  public:
      //! Initialize to default settings.
    Settings();
      //! Copy existing settings
    Settings( const Settings& other );
    ~Settings();
      //! Copy existing settings.
    Settings& operator=( const Settings& other );
  
      //! Set the mapping function for a single element topology.
    void set_mapping_function( const MappingFunction* func );
    void set_mapping_function( const MappingFunction2D* func );
    void set_mapping_function( const MappingFunction3D* func );
      //! Set the mapping function for one or more element topologies.
    void set_mapping_functions( const MappingFunction* const array[], size_t array_size );
    void set_mapping_functions( const MappingFunction2D* const array[], size_t array_size );
    void set_mapping_functions( const MappingFunction3D* const array[], size_t array_size );
      //! Get the mapping function for an element topology.
    const MappingFunction* get_mapping_function( EntityTopology element_type ) const;
    const MappingFunction2D* get_mapping_function_2D( EntityTopology element_type ) const;
    const MappingFunction3D* get_mapping_function_3D( EntityTopology element_type ) const;
    
    enum FixedVertexMode { 
      FIXED_FLAG = 4,   //!< Ask application which vertices are fixed by
                        //!< calling Mesh::vertices_get_fixed_flag
      FIXED_VERTEX = 0, //!< Treat all vertices for which the mesh domain
                        //!< is topologically 0-dimensional as fixed.
      FIXED_CURVE = 1,  //!< Treat all vertices for which the mesh domain
                        //!< is topologically 1-dimensional as fixed.
      FIXED_SURFACE = 2 //!< Treat all vertices for which the corresponding
                        //!< mesh domain has a topological dimension *less than*
                        //!< or equal to 2 as fixed.
    };
      //! Change how Mesquite determines which vertices are fixed.
    void set_fixed_vertex_mode( FixedVertexMode mode );
      //! Get the setting for how Mesquite determines which vertices are fixed.
    FixedVertexMode get_fixed_vertex_mode() const;
    
    enum HigherOrderSlaveMode { 
      SLAVE_NONE,       //!< All higher-order nodes are treated either
                        //!< as fixed or as free variables in the optimization.
 
      SLAVE_ALL,        //!< Treat all non-fixed higher-order nodes as slave 
                        //!< vertices.  The node location will be updated 
                        //!< automatically to the position of the mapping 
                        //!< function evaluated at the logical location of 
                        //!< the node evaluated as if the element did not 
                        //!< contain the node.

      SLAVE_CALCULATED, //!< A utility (e.g. SlaveBoundaryVertices) will be 
                        //!< inserted into the instruction queue that will 
                        //!< determine which vertices are slaved and store 
                        //!< that state as the appropriate bit in the vertex
                        //!< byte.

      SLAVE_FLAG        //!< The results of a call to Mesh::vertices_get_slaved_flag
                        //!< will be used to establish free vs. slaved.
    };
    
      //! Set the slaved higher-order node setting.
    void set_slaved_ho_node_mode( HigherOrderSlaveMode mode );
      //! Get the slaved higher-order node setting.
    HigherOrderSlaveMode get_slaved_ho_node_mode() const;

      /**\brief Generate SIGFPE whenever a floating point exception occurs
       *
       * Generate a FPE signal when overflow, divbyzero, etc. occur
       * during floating-point arithmatic.  This is intended for debugging
       * purposes only, as enabling this will typically result in a 
       * crash when such arithmatic errors occur.
       *
       * If this option is enabled, Mesquite will attempt to set 
       * platform-specific flags such that a SIGFPE is generated for
       * floating point errors while the instruction queue is running.
       * If this option ins disabled, Mesquite will not change the
       * flags.  There is no option to explicitly disable such flags
       * because that is the default behavior on most platforms, and
       * presumably if the application has enabled such flags it makes
       * little sense to disable them while Mesquite is running.
       *
       * This functionality may not be supported on all platforms.  If
       * it is not supported, this option has no effect.
       */
    void trap_floating_point_exception( bool enable );
    bool trap_floating_point_exception() const;

  private:
  
    SettingData* mData;
};


} // namespace Mesquite

#endif
