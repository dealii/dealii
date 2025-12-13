// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test the vectorized approximation of exp function.

#include <deal.II/base/fast_transcendental.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/vectorization.h>

#include <array>
#include <cmath>
#include <iomanip>
#include <type_traits>
#include <typeinfo>

#include "../tests.h"

template <typename T>
struct number_type
{
  using type = T;
};
template <typename N, std::size_t W>
struct number_type<VectorizedArray<N, W>>
{
  using type = N;
};

template <typename Number,
          int                                       pol_degree,
          fast_transcendental::PolynomialEvalScheme pol_eval_scheme>
void
do_specific_test()
{
  using floating_type = typename number_type<Number>::type;

  // The values at indices [0–13] were chosen such that exp(x) evaluates to any
  // finite, non-zero result.
  // To verify correct handling of underflow and overflow, the value at index
  // [14] is chosen so that exp(x) underflows to zero, and the value at index
  // [15] is chosen so that exp(x) overflows to +infinity.
  std::array<floating_type, 16> test_param = {{1.7,
                                               -9.7,
                                               37.0,
                                               5.8,
                                               23.0,
                                               -46.5,
                                               54.9,
                                               -12.2,
                                               42.3,
                                               -78.3,
                                               34.5,
                                               72.0,
                                               -23.4,
                                               34.4,
                                               -812.2,
                                               955.2}};

  deallog << "Polynomial degree " << pol_degree << ". Relative errors: ";

  if constexpr (std::is_floating_point_v<Number>)
    for (unsigned i = 0; i < test_param.size(); ++i)
      {
        if (i > 0)
          deallog << ", ";
        Number v = test_param[i];

        Number exp_approximation = fast_transcendental::
          exp<pol_degree, Number, true, true, pol_eval_scheme>(v);
        Number exp_exact = std::exp(v);
        Number absolute_error =
          v > fast_transcendental::internal::exp_max_abs_x<floating_type> ?
            Number(0.) :
            std::abs(exp_approximation - exp_exact);
        Number relative_error = Number(0.0) == exp_exact ?
                                  absolute_error :
                                  absolute_error / exp_exact;

        relative_error =
          v > fast_transcendental::internal::exp_max_abs_x<floating_type> ?
            exp_approximation :
            relative_error;

        deallog << std::setprecision(10) << relative_error;
      }
  else
    for (unsigned i = 0; i < test_param.size(); i += Number::size())
      {
        if (i > 0)
          deallog << ", ";
        Number v;
        v.load(&test_param[i]);

        Number exp_approximation = fast_transcendental::
          exp<pol_degree, Number, true, true, pol_eval_scheme>(v);
        Number exp_exact = std::exp(v);

        // Computing the relative error for results that are infinity can,
        // depending on the system configuration, trigger a floating-point
        // exception. To avoid this, such values are filtered here and set
        // to 1. Since the test output for this case will later be replaced by
        // the approximated solution anyway, this has no impact on test
        // integrity.
        Number finite_exp_approximation =
          compare_and_apply_mask<SIMDComparison::greater_than>(
            v,
            Number(fast_transcendental::internal::exp_max_abs_x<floating_type>),
            Number(1.0),
            exp_approximation);
        Number finite_exp_exact =
          compare_and_apply_mask<SIMDComparison::greater_than>(
            v,
            Number(fast_transcendental::internal::exp_max_abs_x<floating_type>),
            Number(1.0),
            exp_exact);

        Number absolute_error =
          std::abs(finite_exp_approximation - finite_exp_exact);
        Number relative_error =
          absolute_error / compare_and_apply_mask<SIMDComparison::equal>(
                             exp_exact, Number(0.0), Number(1.0), exp_exact);

        relative_error = compare_and_apply_mask<SIMDComparison::equal>(
          Number(0.0), exp_exact, absolute_error, relative_error);
        relative_error = compare_and_apply_mask<SIMDComparison::greater_than>(
          v,
          Number(fast_transcendental::internal::exp_max_abs_x<floating_type>),
          exp_approximation,
          relative_error);

        deallog << std::setprecision(10) << relative_error;
      }

  deallog << std::endl;
}

template <typename Number, int... ints>
void
do_test(std::integer_sequence<int, ints...>)
{
  using floating_type = typename number_type<Number>::type;

  // The selection of Estrin's scheme for double and Horner's scheme for float
  // is deliberate only to cover both polynomial evaluation mechanisms during
  // testing.
  if constexpr (std::is_same_v<floating_type, double>)
    (do_specific_test<Number,
                      ints,
                      fast_transcendental::PolynomialEvalScheme::estrin>(),
     ...);
  else if constexpr (std::is_same_v<floating_type, float>)
    (do_specific_test<Number,
                      ints,
                      fast_transcendental::PolynomialEvalScheme::horner>(),
     ...);
}

template <typename Number>
void
do_test()
{
  if constexpr (std::is_floating_point_v<Number>)
    deallog << "Type: " << typeid(Number).name() << std::endl;
  else
    deallog << "Type: " << typeid(typename Number::value_type).name()
            << ", width: " << Number::size() << std::endl;
  do_test<Number>(std::integer_sequence<int, 3, 4, 5, 6, 7, 8, 9, 10, 11>{});
}


int
main()
{
  initlog();
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 512
  do_test<dealii::VectorizedArray<float, 16>>();
  do_test<dealii::VectorizedArray<double, 8>>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 256
  do_test<dealii::VectorizedArray<float, 8>>();
  do_test<dealii::VectorizedArray<double, 4>>();
#endif

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS >= 128
  do_test<dealii::VectorizedArray<float, 4>>();
  do_test<dealii::VectorizedArray<double, 2>>();
#endif
  do_test<dealii::VectorizedArray<float, 1>>();
  do_test<dealii::VectorizedArray<double, 1>>();

  do_test<float>();
  do_test<double>();
}
