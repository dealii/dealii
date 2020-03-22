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

    (2011) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file LaplaceWrapper.hpp
 *  \brief Define LaplaceWrapper class
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_LAPLACE_WRAPPER_HPP
#define MSQ_LAPLACE_WRAPPER_HPP

#include "Wrapper.hpp"

namespace MESQUITE_NS {

class LaplaceWrapper : public Wrapper
{
public:
  
  MESQUITE_EXPORT
  LaplaceWrapper();

  /**\brief Specify timeout after which untangler will exit 
   *
   *  Specify a value less than or equal to zero for no limit
   */
  void set_cpu_time_limit( double seconds )
    { maxTime = seconds; }
  double get_cpu_time_limit() const 
    { return maxTime; }

  /**\brief Specify factor by which to minimum distance a vertex must 
   *        move in an iteration to avoid termination of the untangler 
   *
   *  Specify a value less than or equal to zero for no limit.
   *\NOTE Culling cannot be done w/out a limit on vertex movement
   */
  void set_vertex_movement_limit_factor( double f )
    { movementFactor = f; }
  double get_vertex_movement_limit_factor() const
    { return movementFactor; }

  /**\brief Specify maximum number of iterations.  
   * 
   * Specify a value less than or equal to zero for no limit
   */
  void set_iteration_limit( int limit )
    { iterationLimit = limit; }
  int get_iteration_limit() const
    { return iterationLimit; }

  /**\brief Cull vertices based on movement limit */
  inline void enable_culling( bool yesno )
    { doCulling = yesno; }
  inline bool is_culling_enabled() const
    { return doCulling; }
  

  MESQUITE_EXPORT
  ~LaplaceWrapper();

protected:

  MESQUITE_EXPORT
  void run_wrapper( MeshDomainAssoc* mesh_and_domain,
                    ParallelMesh* pmesh,
                    Settings* settings,
                    QualityAssessor* qa,
                    MsqError& err );
  
private:
  
  double maxTime, movementFactor;
  int iterationLimit;
  bool doCulling;
};


} // namespace MESQUITE_NS

#endif
