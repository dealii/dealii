// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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


#ifndef dealii_matrix_free_evaluation_template_factory_internal_h
#define dealii_matrix_free_evaluation_template_factory_internal_h


#include <deal.II/base/config.h>

#ifndef FE_EVAL_FACTORY_DEGREE_MAX
#  define FE_EVAL_FACTORY_DEGREE_MAX 6
#endif

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  template <int degree, typename EvaluatorType, typename... Args>
  bool
  instantiation_helper_run(const unsigned int given_degree,
                           const unsigned int n_q_points_1d,
                           Args &...args)
  {
    if (given_degree == degree)
      {
        if (n_q_points_1d == degree + 1)
          return EvaluatorType::template run<degree, degree + 1>(args...);
        else if (n_q_points_1d == degree + 2)
          return EvaluatorType::template run<degree, degree + 2>(args...);
        else if (n_q_points_1d == degree)
          return EvaluatorType::template run<degree, degree>(args...);
        else if (n_q_points_1d == (3 * degree) / 2 + 1)
          return EvaluatorType::template run<degree, (3 * degree) / 2 + 1>(
            args...);
        else
          // slow path
          return EvaluatorType::template run<-1, 0>(args...);
      }
    else if (degree < FE_EVAL_FACTORY_DEGREE_MAX)
      return instantiation_helper_run<
        (degree < FE_EVAL_FACTORY_DEGREE_MAX ? degree + 1 : degree),
        EvaluatorType>(given_degree, n_q_points_1d, args...);
    else
      // slow path
      return EvaluatorType::template run<-1, 0>(args...);
  }

  template <int degree, typename EvaluatorType, typename... Args>
  bool
  instantiation_helper_degree_run(const unsigned int given_degree,
                                  Args &...args)
  {
    if (given_degree == degree)
      {
        return EvaluatorType::template run<degree, degree + 1>(args...);
      }
    else if (degree < FE_EVAL_FACTORY_DEGREE_MAX)
      return instantiation_helper_degree_run<
        (degree < FE_EVAL_FACTORY_DEGREE_MAX ? degree + 1 : degree),
        EvaluatorType>(given_degree, args...);
    else
      // slow path
      return EvaluatorType::template run<-1, 0>(args...);
  }

} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
