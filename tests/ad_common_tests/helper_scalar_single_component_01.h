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
// Evaluation of a single component (scalar) system using a helper class

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

template <int dim, typename NumberType>
struct FunctionsTestScalar
{
  static NumberType
  psi(const NumberType &s)
  {
    return 4.0 * pow(s, 4);
  }

  static NumberType
  dpsi_ds(const NumberType &s)
  {
    return 16.0 * pow(s, 3);
  }

  static NumberType
  d2psi_ds_ds(const NumberType &s)
  {
    return 48.0 * pow(s, 2);
  }
};

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_scalar()
{
  using ADHelper         = AD::ScalarFunction<dim, ad_type_code, number_t>;
  using ADNumberType     = typename ADHelper::ad_type;
  using ScalarNumberType = typename ADHelper::scalar_type;

  std::cout << "*** Test variables: Scalar dof, "
            << "dim = " << Utilities::to_string(dim) << ", "
            << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  // Values computed from the AD energy function
  ScalarNumberType             psi;
  Vector<ScalarNumberType>     Dpsi;
  FullMatrix<ScalarNumberType> D2psi;

  // Function and its derivatives
  using func_ad = FunctionsTestScalar<dim, ADNumberType>;

  // Setup the variable components and choose a value at which to
  // evaluate the tape
  const FEValuesExtractors::Scalar s_dof(0);
  const unsigned int               n_AD_components = 1;
  ADHelper                         ad_helper(n_AD_components);
  ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                     // default values

  ScalarNumberType s = 1.2;

  // Configure tape
  const int  tape_no = 1;
  const bool is_recording =
    ad_helper.start_recording_operations(tape_no /*material_id*/,
                                         true /*overwrite_tape*/,
                                         true /*keep*/);
  if (is_recording == true)
    {
      ad_helper.register_independent_variable(s, s_dof);

      const ADNumberType s_ad = ad_helper.get_sensitive_variables(s_dof);

      const ADNumberType psi(func_ad::psi(s_ad));

      ad_helper.register_dependent_variable(psi);
      ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);

      std::cout << "Recorded data..." << std::endl;
      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
      std::cout << "s_ad: " << s_ad << std::endl;
      std::cout << "psi: " << psi << std::endl;
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
      s = 1.5;
      ad_helper.activate_recorded_tape(tape_no);
      ad_helper.set_independent_variable(s, s_dof);

      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
    }

  // Compute the function value, gradient and hessian for the new evaluation
  // point
  psi = ad_helper.compute_value();
  ad_helper.compute_gradient(Dpsi);
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    ad_helper.compute_hessian(D2psi);

  // Output the full stored function, gradient vector and hessian matrix
  std::cout << "psi: " << psi << std::endl;
  std::cout << "Dpsi: \n";
  Dpsi.print(std::cout);
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      std::cout << "D2psi: \n";
      D2psi.print_formatted(std::cout, 3, true, 0, "0.0");
    }

  // Extract components of the solution
  const ScalarNumberType dpsi_ds =
    ad_helper.extract_gradient_component(Dpsi, s_dof);

  // Verify the result
  using func = FunctionsTestScalar<dim, ScalarNumberType>;
  static const ScalarNumberType tol =
    1e5 * std::numeric_limits<ScalarNumberType>::epsilon();
  std::cout << "psi:              " << psi << std::endl;
  std::cout << "func::psi(s):     " << func::psi(s) << std::endl;
  Assert(std::abs(psi - func::psi(s)) < tol,
         ExcMessage("No match for function value."));
  std::cout << "dpsi_ds:              " << dpsi_ds << std::endl;
  std::cout << "func::dpsi_ds(s):     " << func::dpsi_ds(s) << std::endl;
  Assert(std::abs(dpsi_ds - func::dpsi_ds(s)) < tol,
         ExcMessage("No match for first derivative."));
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      const ScalarNumberType d2psi_ds_ds =
        ad_helper.extract_hessian_component(D2psi, s_dof, s_dof);
      std::cout << "d2psi_ds_ds:              " << d2psi_ds_ds << std::endl;
      std::cout << "func::d2psi_ds_ds(s):     " << func::d2psi_ds_ds(s)
                << std::endl;
      Assert(std::abs(d2psi_ds_ds - func::d2psi_ds_ds(s)) < tol,
             ExcMessage("No match for second derivative."));
    }

  std::cout << std::endl << std::endl;
}
