// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_matrix_free_evaluation_flags_h
#define dealii_matrix_free_evaluation_flags_h

#include <deal.II/base/config.h>


DEAL_II_NAMESPACE_OPEN



/**
 * @brief The namespace for the EvaluationFlags enum
 *
 * This namespace contains the enum EvaluationFlags used in FEEvaluation
 * to control evaluation and integration of values, gradients, etc..
 */
namespace EvaluationFlags
{
  /**
   * @brief The EvaluationFlags enum
   *
   * This enum contains a set of flags used by FEEvaluation::integrate(),
   * FEEvaluation::evaluate() and others to determine if values, gradients,
   * hessians, or a combination of them is being used.
   */
  enum EvaluationFlags
  {
    /**
     * Do not use or compute anything.
     */
    nothing = 0,
    /**
     * Use or evaluate values.
     */
    values = 0x1,
    /**
     * Use or evaluate gradients.
     */
    gradients = 0x2,
    /**
     * Use or evaluate hessians.
     */
    hessians = 0x4
  };


  /**
   * Global operator which returns an object in which all bits are set which are
   * either set in the first or the second argument. This operator exists since
   * if it did not then the result of the bit-or <tt>operator |</tt> would be an
   * integer which would in turn trigger a compiler warning when we tried to
   * assign it to an object of type UpdateFlags.
   *
   * @ref EvaluationFlags
   */
  inline EvaluationFlags
  operator|(const EvaluationFlags f1, const EvaluationFlags f2)
  {
    return static_cast<EvaluationFlags>(static_cast<unsigned int>(f1) |
                                        static_cast<unsigned int>(f2));
  }



  /**
   * Global operator which sets the bits from the second argument also in the
   * first one.
   *
   * @ref EvaluationFlags
   */
  inline EvaluationFlags &
  operator|=(EvaluationFlags &f1, const EvaluationFlags f2)
  {
    f1 = f1 | f2;
    return f1;
  }


  /**
   * Global operator which returns an object in which all bits are set which are
   * set in the first as well as the second argument. This operator exists since
   * if it did not then the result of the bit-and <tt>operator &</tt> would be
   * an integer which would in turn trigger a compiler warning when we tried to
   * assign it to an object of type UpdateFlags.
   *
   * @ref EvaluationFlags
   */
  inline EvaluationFlags operator&(const EvaluationFlags f1,
                                   const EvaluationFlags f2)
  {
    return static_cast<EvaluationFlags>(static_cast<unsigned int>(f1) &
                                        static_cast<unsigned int>(f2));
  }


  /**
   * Global operator which clears all the bits in the first argument if they are
   * not also set in the second argument.
   *
   * @ref EvaluationFlags
   */
  inline EvaluationFlags &
  operator&=(EvaluationFlags &f1, const EvaluationFlags f2)
  {
    f1 = f1 & f2;
    return f1;
  }

} // namespace EvaluationFlags


DEAL_II_NAMESPACE_CLOSE

#endif
