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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */

/*! \file AveragingQM.hpp
    \brief
Header file for the Mesquite::AveragingQM class

  \author Thomas Leurent
  \author Michael Brewer
  \date   2002-05-01
 */

#ifndef MSQ_AVERAGING_QM_HPP
#define MSQ_AVERAGING_QM_HPP

#include "Mesquite.hpp"
#include "QualityMetric.hpp"

namespace MESQUITE_NS
{
   
     /*! \class AveragingQM
       \brief Averaging functionality for use in quality metrics
     */
   class MsqMeshEntity;
   class PatchData;
   class Vector3D;
   class Matrix3D;
   class MsqError;
   
   class AveragingQM
   {
   public:

     MESQUITE_EXPORT AveragingQM( QualityMetric::AveragingMethod method ) :
       avgMethod(method)
     {}

     virtual ~AveragingQM() {}
     
       /*!Set the averaging method for the quality metric. */
     MESQUITE_EXPORT inline 
     void set_averaging_method( QualityMetric::AveragingMethod method )
      { avgMethod = method; }
      
     MESQUITE_EXPORT inline
     void set_averaging_method( QualityMetric::AveragingMethod method, MsqError& )
      { set_averaging_method(method); }
    
     MESQUITE_EXPORT inline 
     QualityMetric::AveragingMethod get_averaging_method() const
      { return avgMethod; }
     
     //! Return true if the requested averaging scheme is supported
     //! for analytical calculation of gradients.  
     inline bool analytical_average_gradient() 
      { return avgMethod <= QualityMetric::LAST_WITH_GRADIENT; }

     //! Return true if the requested averaging scheme is supported
     //! for analytical calculation of Hessians.  
     inline bool analytical_average_hessian() 
      { return avgMethod <= QualityMetric::LAST_WITH_HESSIAN; }

     //! average_metrics takes an array of length num_values and averages the
     //! contents using averaging method data member avgMethod .
     double average_metrics(const double metric_values[], int num_values,
                            MsqError &err);
                            
     //! Given a list of metric values, calculate the average metric
     //! valude according to the current avgMethod and write into
     //! the passed metric_values array the the value weight/count to
     //! use when averaging gradient vectors for the metric.
     //!\param metric_values : As input, a set of quality metric values
     //!                       to average.  As output, the fraction of
     //!                       the corresponding gradient vector that
     //!                       contributes to the average gradient.
     //!\param num_metric_values The number of values in the passed array.
     double average_metric_and_weights( double metric_values[],
                                        int num_metric_values,
                                        MsqError& err );
     
     /** \brief Average metric values and gradients for per-corner evaluation
      *
      *\param element_type   The element type
      *\param num_corners    The number of corners (e.g. pass 4 for a pyramid
      *                      if the metric couldn't be evaluated for the apex)
      *\param corner_values  An array of metric values, one per element corner
      *\param corner_grads   The corner gradients, 4 for each corner
      *\param vertex_grads   Output.  Gradient at each vertex.
      *\return average metric value for element
      */
     double average_corner_gradients( EntityTopology element_type,
                                  uint32_t fixed_vertices,
                                  unsigned num_corners,
                                  double corner_values[],
                                  const Vector3D corner_grads[],
                                  Vector3D vertex_grads[],
                                  MsqError& err );
     
     /** \brief Average metric values, gradients, and Hessian diagonal 
      *         blocks for per-corner evaluation
      *
      *\param element_type   The element type
      *\param num_corners    The number of corners (e.g. pass 4 for a pyramid
      *                      if the metric couldn't be evaluated for the apex)
      *\param corner_values  An array of metric values, one per element corner
      *\param corner_grads   The corner gradients, 4 for each corner
      *\param corner_hessians The hessians, 10 for each corner
      *\param vertex_grads   Output.  Gradient at each vertex.
      *\param vertex_hessians Output.  Hessian diagonal block for each vertex.
      *\return average metric value for element
      */
      double average_corner_hessian_diagonals( EntityTopology element_type,
                                               uint32_t fixed_vertices,
                                               unsigned num_corners,
                                               const double corner_values[],
                                               const Vector3D corner_grads[],
                                               const Matrix3D corner_hessians[],
                                               Vector3D vertex_grads[],
                                               SymMatrix3D vertex_hessians[],
                                               MsqError& err );
     
     /** \brief Average metric values, gradients, and Hessian diagonal 
      *         blocks for per-corner evaluation
      *
      *\param element_type   The element type
      *\param num_corners    The number of corners (e.g. pass 4 for a pyramid
      *                      if the metric couldn't be evaluated for the apex)
      *\param corner_values  An array of metric values, one per element corner
      *\param corner_grads   The corner gradients, 4 for each corner
      *\param corner_hess_diag The diagonal blocks of the Hessian: 4 for each corner.
      *\param vertex_grads   Output.  Gradient at each vertex.
      *\param vertex_hessians Output.  Hessian diagonal block for each vertex.
      *\return average metric value for element
      */
      double average_corner_hessian_diagonals( EntityTopology element_type,
                                               uint32_t fixed_vertices,
                                               unsigned num_corners,
                                               const double corner_values[],
                                               const Vector3D corner_grads[],
                                               const SymMatrix3D corner_hess_diag[],
                                               Vector3D vertex_grads[],
                                               SymMatrix3D vertex_hessians[],
                                               MsqError& err );
     
     /** \brief Average metric values, gradients, and Hessians for 
      *         per-corner evaluation
      *
      *\param element_type   The element type
      *\param num_corners    The number of corners (e.g. pass 4 for a pyramid
      *                      if the metric couldn't be evaluated for the apex)
      *\param corner_values  An array of metric values, one per element corner
      *\param corner_grads   The corner gradients, 4 for each corner
      *\param corner_hessians The hessians, 10 for each corner
      *\param vertex_grads   Output.  Gradient at each vertex.
      *\param vertex_hessians Output.  Hessians.  Length must be (n*(n+1))/2,
      *                       where n is the number of vertices in the element.
      *\return average metric value for element
      */
      double average_corner_hessians( EntityTopology element_type,
                                     uint32_t fixed_vertices,
                                     unsigned num_corners,
                                     const double corner_values[],
                                     const Vector3D corner_grads[],
                                     const Matrix3D corner_hessians[],
                                     Vector3D vertex_grads[],
                                     Matrix3D vertex_hessians[],
                                     MsqError& err );

  private:
     QualityMetric::AveragingMethod avgMethod;
   };
   
} //namespace


#endif // MSQ_AVERAGING_QM_HPP
