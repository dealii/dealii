/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
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
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file TMPQualityMetric.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_TMP_QUALITY_METRIC_HPP
#define MSQ_TMP_QUALITY_METRIC_HPP

#include "Mesquite.hpp"
#include "ElemSampleQM.hpp"
#include "MsqMatrix.hpp"

namespace MESQUITE_NS {

class TargetCalculator;
class WeightCalculator;
class TRel2DMetric;
class TRel3DMetric;
class NodeSet;
class Mesh;
class MeshDomain;
class Settings;

/**\brief Compare targets to mapping function Jacobian matrices
 *
 * Base class for various TMP QualityMetric implementations
 */
class TMPQualityMetric : public ElemSampleQM
{
public:

  /**
   *\param tc   The target calculator 
   *\param wc   The weight calculator
   */
  TMPQualityMetric( TargetCalculator* tc,
                    WeightCalculator* wc ) 
    : targetCalc(tc),
      weightCalc(wc)
   {}
  
  MESQUITE_EXPORT virtual 
  int get_negate_flag() const;
  
  MESQUITE_EXPORT virtual
  void get_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles, 
                        bool free_vertices_only,
                        MsqError& err );
  
  MESQUITE_EXPORT static
  void get_patch_evaluations( PatchData& pd, 
                        std::vector<size_t>& handles, 
                        bool free_vertices_only,
                        MsqError& err );
  
  MESQUITE_EXPORT virtual 
  void get_element_evaluations( PatchData& pd, size_t elem_index,
                                std::vector<size_t>& handles,
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
  
  MESQUITE_EXPORT
  void set_target_calculator( TargetCalculator* tc ) { targetCalc = tc; }
  MESQUITE_EXPORT
  void set_weight_calculator( WeightCalculator* wc ) { weightCalc = wc; }
  MESQUITE_EXPORT
  TargetCalculator* get_target_calculator() const { return targetCalc; }
  MESQUITE_EXPORT
  WeightCalculator* get_weight_calculator() const { return weightCalc; }
    
  MESQUITE_EXPORT
  virtual void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                                 const Settings* settings,
                                 MsqError& err );
  
protected:
 
  MESQUITE_EXPORT virtual
  bool evaluate_internal( PatchData& pd,
                 size_t handle,
                 double& value,
                 size_t* indices,
                 size_t& num_indices,
                 MsqError& err ) = 0;

  MESQUITE_EXPORT
  bool evaluate_surface_common( // input:
                                PatchData& pd,
                                Sample sample,
                                size_t element_index,
                                const NodeSet& bits,
                                // output:
                                size_t* indices, 
                                size_t& num_indices,
                                MsqVector<2>* derivs,
                                MsqMatrix<2,2>& W,
                                MsqMatrix<2,2>& A,
                                MsqMatrix<3,2>& S_a_transpose_Theta,
                                MsqError& err );
                                

   MESQUITE_EXPORT
 void weight( PatchData& pd,
               Sample sample,
               size_t elem,
               int num_points,
               double& value,
               Vector3D* grad,
               SymMatrix3D* diag,
               Matrix3D* hess,
               MsqError& err );
  
  enum { MAX_ELEM_NODES = 27 };
  size_t mIndices[MAX_ELEM_NODES];
  std::vector< MsqMatrix<2,2> > hess2d;
  MsqVector<3> mDerivs3D[MAX_ELEM_NODES];
  MsqVector<2> mDerivs2D[MAX_ELEM_NODES];

  TargetCalculator* targetCalc;
  
private:
  WeightCalculator* weightCalc;
};

} // namespace Mesquite

#endif
