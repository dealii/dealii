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


/** \file UntangleWrapper.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_UNTANGLE_WRAPPER_HPP
#define MSQ_UNTANGLE_WRAPPER_HPP

#include "Mesquite.hpp"
#include "Wrapper.hpp"

namespace MESQUITE_NS {

/**\brief Wrapper that implements several TMP-based untanglers */
class UntangleWrapper : public Wrapper {
  
public:

  /**\brief Which quality metric to use */
  enum UntangleMetric { 
    BETA, //!< Use TMP UntangleBeta metric
    SIZE, //!< Use UntangleMu(TargetSize}
    SHAPESIZE //!< Use UntangleMu(TargetShapeSize}
  };

  /**\brief Specify which untangle metric to use */
  MESQUITE_EXPORT
  void set_untangle_metric( UntangleMetric metric );

  /**\brief Specify constant value for untangle metric */
  MESQUITE_EXPORT
  void set_metric_constant( double value );

  /**\brief Specify timeout after which untangler will exit */
  MESQUITE_EXPORT
  void set_cpu_time_limit( double seconds );

  /**\brief Specify max number of outer iterations after which untangler will exit */
  MESQUITE_EXPORT
  void set_outer_iteration_limit( int maxIt );

  /**\brief Specify factor by which to minimum distance a vertex must 
   *        move in an iteration to avoid termination of the untangler */
  MESQUITE_EXPORT
  void set_vertex_movement_limit_factor( double f );

  MESQUITE_EXPORT
  UntangleWrapper();

  MESQUITE_EXPORT
  UntangleWrapper(UntangleMetric m);

  MESQUITE_EXPORT
  ~UntangleWrapper();

  /**\brief Check if vertex culling will be used */
  inline bool is_culling_enabled() const
    { return doCulling; }
  
  /**\brief Enable vertex culling */
  inline void enable_culling( bool yesno )
    { doCulling = yesno; }

  /**\brief Check if a Jacobi optimization strategy will be used */
  inline bool is_jacobi_optimization() const
    { return doJacobi; }
  
  /**\brief Check if a Gauss optimization strategy will be used */
  inline bool is_gauss_optimization() const
    { return !doJacobi; }
    
  /**\brief Use a Jacobi optimization strategy */
  inline void do_jacobi_optimization()
    { doJacobi = true; }
  
  /**\brief Use a Gauss optimization strategy */
  inline void do_gauss_optimization()
    { doJacobi = false; }

protected:

  MESQUITE_EXPORT
  void run_wrapper( MeshDomainAssoc* mesh_and_domain,
                    ParallelMesh* pmesh,
                    Settings* settings,
                    QualityAssessor* qa,
                    MsqError& err );

private:

  UntangleMetric qualityMetric;
  double maxTime, movementFactor, metricConstant;
  int maxIterations;
  bool doCulling, doJacobi;
};


} // namespace MESQUITE_NS

#endif
