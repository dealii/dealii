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


/** \file ElemSampleQM.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_ELEM_SAMPLE_QM_HPP
#define MSQ_ELEM_SAMPLE_QM_HPP

#include "Mesquite.hpp"
#include "QualityMetric.hpp"
#include "Sample.hpp"

namespace MESQUITE_NS {

/**\brief Base type for metrics evaluated at several 
 *        sample points within each element.
 *
 * This class defines an interface for metrics that are evaluated
 * at multiple sample points within each element.  This interface is
 * used by ElemAverageQM and is the base for ElemCornerQM and
 * ElemTargetQM.
 */
class ElemSampleQM : public QualityMetric
{
public:
  
  MESQUITE_EXPORT virtual ~ElemSampleQM();
  
  MESQUITE_EXPORT virtual MetricType get_metric_type() const
    { return ELEMENT_BASED; }
  
    /**\brief Get evaluation point handles for a given element
     *
     * Similar to QualityMetric::get_evaluations, this method returns
     * a list of handles corresponding to sample points at which the
     * metric may be evaluated.  While QualityMetric::get_evaluations
     * returns sample points for all elements in a PatchData, this
     * method returns only the subset corresponding to a single element.
     */
  MESQUITE_EXPORT virtual 
  void get_element_evaluations( PatchData& pd, size_t elem_index,
                                std::vector<size_t>& handles,
                                MsqError& err ) = 0;

  /** Misc constants used in defining how element index, side dimension,
   *  and side number are packed into a single handle describing a logical
   *  location in the patch at which a sample-based metric is to be 
   *  evaluated.
   */
#ifndef _MSC_VER
  enum {
    /** the number of bits in a handle that are used to store element index */
    ELEM_INDEX_BITS = sizeof(size_t)*8 - Sample::NUMBER_PACKED_BITS,
    /** the maximum number of elements in a PatchData without overflowing handle space */
    MAX_ELEM_PER_PATCH = ((size_t)1)<<ELEM_INDEX_BITS,
    /** Mask to remove sample bits from handle */
    ELEM_SAMPLE_MASK = MAX_ELEM_PER_PATCH - 1
  };
#else /* MS Visual C compiler broken for 64-bit enums */
  static const size_t ELEM_INDEX_BITS = sizeof(size_t)*8 - Sample::NUMBER_PACKED_BITS;
  static const size_t MAX_ELEM_PER_PATCH = ((size_t)1)<<ELEM_INDEX_BITS;
  static const size_t ELEM_SAMPLE_MASK = MAX_ELEM_PER_PATCH - 1;
#endif
  
  inline static size_t handle( Sample sample, size_t index )
    { return (sample.pack() << ELEM_INDEX_BITS) | index; }
  
  inline static Sample sample( size_t handle ) 
    { return Sample(handle >> ELEM_INDEX_BITS); }
  
  inline static size_t elem( size_t handle ) 
    { return handle & ELEM_SAMPLE_MASK; }
};

} // namespace Mesquite

#endif
