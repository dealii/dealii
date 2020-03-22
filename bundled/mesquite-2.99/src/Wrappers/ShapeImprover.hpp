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


/** \file ShapeImprover.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_SHAPE_IMPROVER_HPP
#define MSQ_SHAPE_IMPROVER_HPP

#include "Mesquite.hpp"
#include "Wrapper.hpp"

namespace MESQUITE_NS {

  /**\brief Wrapper that implements TMP-based shape improvement
     */
  class ShapeImprover : public Wrapper {
     
  public:  

    MESQUITE_EXPORT
    ShapeImprover();

    /*\brief Set limit on seconds of CPU time
     *
     * Default is no limit on CPU time.
     */
    MESQUITE_EXPORT
    void set_cpu_time_limit( double seconds ); 

    /*\brief Set mean edge length factor for termination
     *
     * Optimization will cease when, for a given iteration,
     * no vertex is moved more than beta * (mean - sigma)
     * where mean is the average edge length of the initial mesh
     * and sigma is the standard deviation of the edge lengths
     * in the initial mesh.  beta is the value passed to this
     * function.  beta must be greater than zero and less than one.
     * The default value for beta is 10 percent of the meshes 
     * minimum edge length.
     */
    MESQUITE_EXPORT
    void set_vertex_movement_limit_factor( double beta );
    
    MESQUITE_EXPORT
    void set_parallel_iterations( int count );
    

  protected:
  
    MESQUITE_EXPORT
    void run_wrapper( MeshDomainAssoc* mesh_and_domain,
                      ParallelMesh* pmesh,
                      Settings* settings,
                      QualityAssessor* qa,
                      MsqError& err );
    
      
  private:

    double maxTime, mBeta;
    int parallelIterations;
  };
  
  
} // namespace MESQUITE_NS

#endif
