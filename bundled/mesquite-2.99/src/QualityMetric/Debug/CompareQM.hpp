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


/** \file CompareQM.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_COMPARE_QM_HPP
#define MSQ_COMPARE_QM_HPP

#include "Mesquite.hpp"
#include "QualityMetric.hpp"
#include "SimpleStats.hpp"

namespace MESQUITE_NS {

/**\brief Compare values for two supposedly equivalent quality metrics
 *
 * Evaluate two different quality metrics during the evaluation,
 * comparing the results with each other and passing the results
 * of the primary one on to the objective function.
 *
 * Both metrics must be of the same "type", meaning that they must be
 * evaluated at the same sample locations in the mesh.  For example,
 * an error will be generated if one of the metrics is vertex based
 * and one is element based.  Further, both metrics must use the 
 * same handle values to indicate their list of sample locations 
 * such that the evaluation for a given same handle dependes on the
 * same vertices for both metrics.
 */
class CompareQM : public QualityMetric
{
public:
  MESQUITE_EXPORT
  CompareQM( QualityMetric* primary,
             QualityMetric* other,
             const char* primary_name = 0,
             const char* other_name = 0 );

  MESQUITE_EXPORT
  void abort_on_mismatch( double tolerance_factor = 1e-6 );

  MESQUITE_EXPORT
  void do_not_abort();

  MESQUITE_EXPORT
  bool will_abort_on_mismatch() const;
  
  MESQUITE_EXPORT
  void print_stats() const;
  
  MESQUITE_EXPORT virtual
  ~CompareQM();
     
  MESQUITE_EXPORT virtual 
  MetricType get_metric_type() const;
  
  MESQUITE_EXPORT virtual
  std::string get_name() const;
  
  MESQUITE_EXPORT virtual
  int get_negate_flag() const;
  
  MESQUITE_EXPORT virtual
  void get_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles, 
                        bool free_vertices_only,
                        MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate( PatchData& pd, 
                 size_t handle, 
                 double& value, 
                 MsqError& err );


  MESQUITE_EXPORT virtual
  bool evaluate_with_indices( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_gradient( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_Hessian_diagonal( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 std::vector<SymMatrix3D>& Hessian_diagonal,
                 MsqError& err );

  MESQUITE_EXPORT virtual
  bool evaluate_with_Hessian( PatchData& pd,
                 size_t handle,
                 double& value,
                 std::vector<size_t>& indices,
                 std::vector<Vector3D>& gradient,
                 std::vector<Matrix3D>& Hessian,
                 MsqError& err );
private:
  
  double epsilon( double val1, double val2 );
  
  bool check_valid( size_t handle, bool valid1, bool valid2 );


    /** Compare values, add to stats, etc. */
  void check_value( size_t handle, double value1, double value2 );
  
    /** Check that two index lists are equivalent and return
     *  in \c map_out for each index in \c idx1 the position
     *  of the same index in \c idx2.
     */
  void check_indices( size_t handle, 
                      const std::vector<size_t>& idx1,
                      const std::vector<size_t>& idx2, 
                      std::vector<size_t>& map_out,
                      MsqError& err );
                      
  void check_grad( size_t handle, 
                   const std::vector<size_t>& indices,
                   const std::vector<size_t>& index_map,
                   const std::vector<Vector3D>& grad1,
                   const std::vector<Vector3D>& grad2 );
                   
  void check_hess_diag( size_t handle, 
                        const std::vector<size_t>& indices,
                        const std::vector<size_t>& index_map,
                        const std::vector<SymMatrix3D>& hess1,
                        const std::vector<SymMatrix3D>& hess2 );
                   
  void check_hess( size_t handle, 
                   const std::vector<size_t>& indices,
                   const std::vector<size_t>& index_map,
                   const std::vector<Matrix3D>& hess1,
                   const std::vector<Matrix3D>& hess2 );

  void index_mismatch( size_t handle,
                       const std::vector<size_t>& idx1,
                       const std::vector<size_t>& idx2,
                       MsqError& err );

  struct GradStat {
    SimpleStats x, y, z;
    void add( Vector3D grad );
    void add_diff( Vector3D grad1, Vector3D grad2 );
  };
  struct HessStat {
    SimpleStats xx, xy, xz, yy, yz, zz;
    void add_diag( Matrix3D hess );
    void add_diag( SymMatrix3D hess );
    void add_diag_diff( Matrix3D hess1, Matrix3D hess2 );
    void add_diag_diff( SymMatrix3D hess1, SymMatrix3D hess2 );
    void add_nondiag( Matrix3D hess );
    void add_nondiag_diff( Matrix3D hess1, Matrix3D hess2 );
  };

  std::string primaryName, otherName;
  QualityMetric *primaryMetric, *otherMetric;
  bool abortOnMismatch;
  double toleranceFactor;
  SimpleStats valPrimary, valOther, valDiff;
  GradStat gradPrimary, gradOther, gradDiff;
  HessStat hessPrimary, hessOther, hessDiff;

};

} // namespace MESQUITE_NS

#endif
