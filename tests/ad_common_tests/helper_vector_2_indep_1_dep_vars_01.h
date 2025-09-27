// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Header file:
// Evaluation of a vector of 1 dependent and 2 independent variables
// using a helper class

#include <deal.II/base/logstream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

namespace AD = dealii::Differentiation::AD;

// Function and its derivatives
template <int dim, typename NumberType>
struct FunctionsTestSquare
{
  static NumberType
  f0(const NumberType &s0, const NumberType &s1)
  {
    return 2.0 * pow(s0, 4) * pow(s1, 3);
  };

  static NumberType
  df0_ds0(const NumberType &s0, const NumberType &s1)
  {
    return 8.0 * pow(s0, 3) * pow(s1, 3);
  };

  static NumberType
  df0_ds1(const NumberType &s0, const NumberType &s1)
  {
    return 6.0 * pow(s0, 4) * pow(s1, 2);
  };
};

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_AD_vector_jacobian()
{
  using ADHelper         = AD::VectorFunction<dim, ad_type_code, number_t>;
  using ADNumberType     = typename ADHelper::ad_type;
  using ScalarNumberType = typename ADHelper::scalar_type;

  std::cout << "*** Test variables: Variables: 2 independent, 1 dependent, "
            << "dim = " << Utilities::to_string(dim) << ", "
            << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  // Function and its derivatives
  using func_ad = FunctionsTestSquare<dim, ADNumberType>;

  // Setup the variable components and choose a value at which to
  // evaluate the tape
  const unsigned int n_vars_indep = 2;
  const unsigned int n_vars_dep   = 1;
  ADHelper           ad_helper(n_vars_indep, n_vars_dep);
  ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                     // default values

  std::vector<ScalarNumberType> s(n_vars_indep);
  s[0] = 3.1;
  s[1] = 5.9;

  // Configure tape
  const int  tape_no = 1;
  const bool is_recording =
    ad_helper.start_recording_operations(tape_no /*material_id*/,
                                         true /*overwrite_tape*/,
                                         true /*keep*/);
  if (is_recording == true)
    {
      ad_helper.register_independent_variables(s);

      const std::vector<ADNumberType> s_ad =
        ad_helper.get_sensitive_variables();

      std::vector<ADNumberType> f_ad(n_vars_dep, ADNumberType(0.0));
      f_ad[0] = func_ad::f0(s_ad[0], s_ad[1]);

      ad_helper.register_dependent_variables(f_ad);
      ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);

      std::cout << "Taped data..." << std::endl;
      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
      std::cout << "s_ad: ";
      for (unsigned int i = 0; i < n_vars_indep; ++i)
        std::cout << s_ad[i] << (i < (n_vars_indep - 1) ? "," : "");
      std::cout << std::endl;
      std::cout << "f_ad: ";
      for (unsigned int i = 0; i < n_vars_dep; ++i)
        std::cout << f_ad[i] << (i < (n_vars_dep - 1) ? "," : "");
      std::cout << std::endl;
    }
  else
    {
      Assert(is_recording == true, ExcInternalError());
    }

  // Do some work :-)
  // Set a new evaluation point
  if (AD::ADNumberTraits<ADNumberType>::is_taped == true)
    {
      std::cout
        << "Using tape with different values for independent variables..."
        << std::endl;
      s[0] = 4.9;
      s[1] = 0.87;
      ad_helper.activate_recorded_tape(tape_no);
      ad_helper.set_independent_variables(s);

      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
    }

  // Compute the function values and their jacobian for the new evaluation point
  Vector<ScalarNumberType>     funcs(n_vars_dep);
  FullMatrix<ScalarNumberType> Dfuncs(n_vars_dep, n_vars_indep);
  ad_helper.compute_values(funcs);
  ad_helper.compute_jacobian(Dfuncs);

  // Output the full stored function, gradient vector and hessian matrix
  std::cout << "funcs: \n";
  funcs.print(std::cout);
  std::cout << "Dfuncs: \n";
  Dfuncs.print_formatted(std::cout, 3, true, 0, "0.0");

  // Verify the result
  using func = FunctionsTestSquare<dim, ScalarNumberType>;
  static const ScalarNumberType tol =
    1e5 * std::numeric_limits<ScalarNumberType>::epsilon();
  std::cout << "funcs[0]: " << funcs[0]
            << "\t func::f0(s[0],s[1])): " << func::f0(s[0], s[1]) << std::endl;
  std::cout << "Dfuncs[0][1]: " << Dfuncs[0][0]
            << "\t func::df0_ds0(s[0],s[1])): " << func::df0_ds0(s[0], s[1])
            << std::endl;
  std::cout << "Dfuncs[0][1]: " << Dfuncs[0][1]
            << "\t func::df0_ds1(s[0],s[1])): " << func::df0_ds1(s[0], s[1])
            << std::endl;
  Assert(std::abs(funcs[0] - func::f0(s[0], s[1])) < tol,
         ExcMessage("No match for function 1 value."));
  Assert(std::abs(Dfuncs[0][0] - func::df0_ds0(s[0], s[1])) < tol,
         ExcMessage("No match for function 1 first derivative.."));
  Assert(std::abs(Dfuncs[0][1] - func::df0_ds1(s[0], s[1])) < tol,
         ExcMessage("No match for function 1 first derivative.."));
}
