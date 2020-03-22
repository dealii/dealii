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

    (2006) kraftche@cae.wisc.edu    
   
  ***************************************************************** */

/*! \file QualityMetric.hpp
    \brief
Header file for the Mesquite::QualityMetric class

  \author Thomas Leurent
  \author Michael Brewer
  \date   2002-05-01
 */

#ifndef QualityMetric_hpp
#define QualityMetric_hpp

#include <cmath>
#include <vector>
#include <algorithm>

#include "Mesquite.hpp"
#include "Vector3D.hpp"
#include "Matrix3D.hpp"

#ifdef _MSC_VER
   typedef unsigned uint32_t;
#elif defined(MSQ_HAVE_STDINT_H)
#  include <stdint.h>
#elif defined(MSQ_HAVE_INTTYPES_H)
#  include <inttypes.h>
#endif

namespace MESQUITE_NS
{
   
     /*! \class QualityMetric
       \brief Base class for concrete quality metrics.
     */
   class PatchData;
   class MsqMeshEntity;
   class Mesh;
   class MeshDomain;
   class MeshDomainAssoc;
   class Settings;
   
   class QualityMetric
   {
   protected:

     QualityMetric( ) : 
      keepFiniteDiffEps(false), 
      haveFiniteDiffEps(false) 
      {}

   public:

     enum MetricType
     {
        VERTEX_BASED,  /**< Iterate over vertices to evaluate metric. */
        ELEMENT_BASED  /**< Iterate over elements to evaluate metric. */
     };

     MESQUITE_EXPORT virtual ~QualityMetric()
      {}
     
     MESQUITE_EXPORT virtual MetricType get_metric_type() const = 0;
     
     MESQUITE_EXPORT virtual std::string get_name() const = 0;

      //! 1 if metric should be minimized, -1 if metric should be maximized.
     MESQUITE_EXPORT virtual int get_negate_flag() const = 0;
     
      /**\brief Get locations at which metric can be evaluated
       *
       * Different metrics are evaluated for different things within
       * a patch.  For example, an element-based metric will be evaluated
       * once for each element in patch, a vertex-based metric once for 
       * each veretx, and a target/sample-point based metric will be 
       * evaluated once for each samle point in each element.  This method
       * returns a list of handles, one for each location in the patch
       * at which the metric can be evaluated.  The handle values are used
       * as input to the evaluate methods.
       *\param pd       The patch
       *\param handles  Output list of handles
       *\param free_vertices_only If true, only pass back evaluation points
       *         that depend on at least one free vertex.
       */
     MESQUITE_EXPORT virtual
     void get_evaluations( PatchData& pd, 
                           std::vector<size_t>& handles, 
                           bool free_vertices_only,
                           MsqError& err ) = 0;
     
     
      /**\brief Get locations at which metric can be evaluated for
       *        use in BCD intialization and QualityAssessor.
       *
       * For element-based, sample-based, and vertex-based metrics,
       * this function is the same as get_evaluations.  For edge-based
       * metrics it returns only a subset of the results for get_evaluations
       * such that each edge in the mesh is visited only once even though
       * it would normally be visited twice when iterating over patches
       * of the mesh.  This assumes that no vertex occurs in more than one
       * patch without its MSQ_PATCH_FIXED flag set.  This assumption is true for
       * both element-on-vertex and global patches.
       *\param pd       The patch
       *\param handles  Output list of handles
       *\param free_vertices_only If true, only pass back evaluation points
       *         that depend on at least one free vertex.
       */
     MESQUITE_EXPORT virtual
     void get_single_pass( PatchData& pd,
                           std::vector<size_t>& handles, 
                           bool free_vertices_only,
                           MsqError& err );
     
     /**\brief Get metric value at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd     The patch.
      *\param handle The location in the patch (as passed back from get_evaluations).
      *\param value  The output metric value.
      */
     MESQUITE_EXPORT virtual
     bool evaluate( PatchData& pd, 
                    size_t handle, 
                    double& value, 
                    MsqError& err ) = 0;
     
     /**\brief Get metric value at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd      The patch.
      *\param handle  The location in the patch (as passed back from get_evaluations).
      *\param value   The output metric value.
      *\param indices The free vertices that the evaluation is a function
      *               of, specified as vertex indices in the PatchData.
      */
     MESQUITE_EXPORT virtual
     bool evaluate_with_indices( PatchData& pd,
                    size_t handle,
                    double& value,
                    std::vector<size_t>& indices,
                    MsqError& err ) = 0;
     
     /**\brief Get metric value and gradient at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd      The patch.
      *\param handle  The location in the patch (as passed back from get_evaluations).
      *\param value   The output metric value.
      *\param indices The free vertices that the evaluation is a function
      *               of, specified as vertex indices in the PatchData.
      *\param gradient The gradient of the metric as a function of the
      *               coordinates of the free vertices passed back in
      *               the indices list.
      */
     MESQUITE_EXPORT virtual
     bool evaluate_with_gradient( PatchData& pd,
                    size_t handle,
                    double& value,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& gradient,
                    MsqError& err );
     
     /**\brief Get metric value and gradient at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd      The patch.
      *\param handle  The location in the patch (as passed back from get_evaluations).
      *\param value   The output metric value.
      *\param indices The free vertices that the evaluation is a function
      *               of, specified as vertex indices in the PatchData.
      *\param gradient The gradient of the metric as a function of the
      *               coordinates of the free vertices passed back in
      *               the indices list.
      *\param Hessian_diagonal The 3x3 blocks along the diagonal of
      *               the Hessian matrix.
      */
     MESQUITE_EXPORT virtual
     bool evaluate_with_Hessian_diagonal( PatchData& pd,
                    size_t handle,
                    double& value,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& gradient,
                    std::vector<SymMatrix3D>& Hessian_diagonal,
                    MsqError& err );
     
     /**\brief Get metric value and deravitives at a logical location in the patch.
      *
      * Evaluate the metric at one location in the PatchData.
      *\param pd      The patch.
      *\param handle  The location in the patch (as passed back from get_evaluations).
      *\param value   The output metric value.
      *\param indices The free vertices that the evaluation is a function
      *               of, specified as vertex indices in the PatchData.
      *\param gradient The gradient of the metric as a function of the
      *               coordinates of the free vertices passed back in
      *               the indices list.
      *\param Hessian The Hessian of the metric as a function of the 
      *               coordinates. The Hessian is passed back as the
      *               upper-triangular portion of the matrix in row-major
      *               order, where each Matrix3D is the portion of the
      *               Hessian with respect to the vertices at the
      *               corresponding positions in the indices list.
      */
     MESQUITE_EXPORT virtual
     bool evaluate_with_Hessian( PatchData& pd,
                    size_t handle,
                    double& value,
                    std::vector<size_t>& indices,
                    std::vector<Vector3D>& gradient,
                    std::vector<Matrix3D>& Hessian,
                    MsqError& err );

       //!Escobar Barrier Function for Shape and Other Metrics
       // det = signed determinant of Jacobian Matrix at a Vertex
       // delta = scaling parameter
     static inline double vertex_barrier_function(double det, double delta) 
            { return 0.5*(det+sqrt(det*det+4*delta*delta)); }
  //protected:

      /** \brief Remove from vector any gradient terms corresponding 
       *         to a fixed vertex.
       *
       * Remove terms from vector that correspond to fixed vertices.
       *\param type            Element type
       *\param fixed_vertices  Bit flags, one per vertex, 1 if
       *                       vertex is fixed.
       *\param gradients       Array of gradients
       */
     MESQUITE_EXPORT  
	 static void remove_fixed_gradients( EntityTopology type, 
                                          uint32_t fixed_vertices, 
                                          std::vector<Vector3D>& gradients );

      /** \brief Remove from vectors any gradient terms and hessian
       *         diagonal blcoks corresponding to a fixed vertex.
       *
       * Remove terms from vector that correspond to fixed vertices.
       *\param type            Element type
       *\param fixed_vertices  Bit flags, one per vertex, 1 if
       *                       vertex is fixed.
       *\param gradients       Array of gradients
       *\param hess_diagonal_blocks   Array of diagonal blocks of Hessian matrix.
       */
     MESQUITE_EXPORT  
	 static void remove_fixed_diagonals( EntityTopology type, 
                                          uint32_t fixed_vertices, 
                                          std::vector<Vector3D>& gradients,
                                          std::vector<SymMatrix3D>& hess_diagonal_blocks );

      /** \brief Remove from vector any Hessian blocks corresponding 
       *         to a fixed vertex.
       *
       * Remove blocks from vector that correspond to fixed vertices.
       *\param type            Element type
       *\param fixed_vertices  Bit flags, one per vertex, 1 if
       *                       vertex is fixed.
       *\param hessians        Array of Hessian blocks (upper trianguler, row-major)
       */
     MESQUITE_EXPORT  
	 static void remove_fixed_hessians ( EntityTopology type, 
                                          uint32_t fixed_vertices, 
                                          std::vector<Matrix3D>& hessians );
     
     /** \brief Convert fixed vertex format from list to bit flags
      *
      * Given list of pointers to fixed vertices as passed to
      * evaluation functions, convert to bit flag format used
      * for many utility functions in this class.  Bits correspond
      * to vertices in the canonical vertex order, beginning with
      * the least-significant bit.  The bit is cleared for free
      * vertices and set (1) for fixed vertices.
      */
     MESQUITE_EXPORT  
	 static uint32_t fixed_vertex_bitmap( PatchData& pd, 
                                           const MsqMeshEntity* elem,
                                           std::vector<size_t>& free_indices );
      
     
     //! takes an array of coefficients and an array of metrics (both of length num_value)
     //! and averages the contents using averaging method 'method'.
     MESQUITE_EXPORT 
	 double weighted_average_metrics(const double coef[],
                                    const double metric_values[],
                                    const int& num_values, MsqError &err);

       /*!AveragingMethod allows you to set how the quality metric values
         attained at each sample point will be averaged together to produce
         a single metric value for an element.
       */
     enum AveragingMethod
     {
        LINEAR,                 //!< the linear average
        RMS,                    //!< the root-mean-squared average
        HMS,                    //!< the harmonic-mean-squared average
        SUM,                    //!< the sum of the values
        SUM_SQUARED,            //!< the sum of the squares of the values
        HARMONIC,               //!< the harmonic average
        LAST_WITH_HESSIAN=HARMONIC,
        MINIMUM,                //!< the minimum value
        MAXIMUM,                //!< the maximum value
        GEOMETRIC,              //!< the geometric average
        LAST_WITH_GRADIENT=GEOMETRIC,
        STANDARD_DEVIATION,     //!< the standard deviation squared of the values
        MAX_OVER_MIN,           //!< the maximum value minus the minum value
        MAX_MINUS_MIN,          //!< the maximum value divided by the minimum value
        SUM_OF_RATIOS_SQUARED   //!< (1/(N^2))*(SUM (SUM (v_i/v_j)^2))
     };
     
      //!\brief Called at start of instruction queue processing
      //!
      //! Do any preliminary global initialization, consistency checking,
      //! etc.  Default implementation does nothing.
     MESQUITE_EXPORT virtual 
     void initialize_queue( MeshDomainAssoc* mesh_and_domain,
                            const Settings* settings,
                            MsqError& err );

  private:
     int feasible;
     
     std::vector<Matrix3D> tmpHess;
     bool keepFiniteDiffEps; //!< True if gradient finite difference
                             //!< calculation should set \c finiteDiffEps
     bool haveFiniteDiffEps; //!< True if finite difference Hessian code
                             //!< has calculated \c finiteDiffEps 
     double finiteDiffEps; //!< Location for finite difference Hessian code
                           //!< to store this value so that it doesn't need
                           //!< to be recalculated if the gradient calculation
                           //!< is also finite difference
   };


} //namespace


#endif // QualityMetric_hpp
