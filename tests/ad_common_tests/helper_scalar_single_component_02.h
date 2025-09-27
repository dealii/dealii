// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2021 by the deal.II authors
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
// Evaluation of a single component (vector) system using a helper class

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
struct FunctionsTestVector
{
  static NumberType
  psi(const Tensor<1, dim, NumberType> &v)
  {
    // Potential memory corruption in adtl::adouble...
    // Probably need to zero some temporary value
    // in Tensor::operator* or Tensor::contract()
    return 2.0 * (v * v);
    // return 2.0*contract(v,v);
    // return 2.0*scalar_product(v,v);
  }

  static Tensor<1, dim, NumberType>
  dpsi_dv(const Tensor<1, dim, NumberType> &v)
  {
    return 4.0 * v;
  }

  static Tensor<2, dim, NumberType>
  d2psi_dv_dv(const Tensor<1, dim, NumberType> &v)
  {
    const SymmetricTensor<2, dim, NumberType> I(unit_symmetric_tensor<dim>());
    const Tensor<2, dim, NumberType>          I_ns(I);
    return 4.0 * I_ns;
  }
};

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_vector()
{
  using ADHelper         = AD::ScalarFunction<dim, ad_type_code, number_t>;
  using ADNumberType     = typename ADHelper::ad_type;
  using ScalarNumberType = typename ADHelper::scalar_type;

  std::cout << "*** Test variables: Vector dof, "
            << "dim = " << Utilities::to_string(dim) << ", "
            << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  // Values computed from the AD energy function
  ScalarNumberType             psi;
  Vector<ScalarNumberType>     Dpsi;
  FullMatrix<ScalarNumberType> D2psi;

  // Function and its derivatives
  using func_ad = FunctionsTestVector<dim, ADNumberType>;

  // Setup the variable components and choose a value at which to
  // evaluate the tape
  const FEValuesExtractors::Vector v_dof(0);
  const unsigned int               n_AD_components = dim;
  ADHelper                         ad_helper(n_AD_components);
  ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                     // default values

  Tensor<1, dim, ScalarNumberType> v;
  for (unsigned int i = 0; i < dim; ++i)
    v[i] = 1.0 + i;

  const int  tape_no = 1;
  const bool is_recording =
    ad_helper.start_recording_operations(tape_no /*material_id*/,
                                         true /*overwrite_tape*/,
                                         true /*keep*/);
  if (is_recording == true)
    {
      ad_helper.register_independent_variable(v, v_dof);

      const Tensor<1, dim, ADNumberType> v_ad =
        ad_helper.get_sensitive_variables(v_dof);

      const ADNumberType psi = func_ad::psi(v_ad);

      ad_helper.register_dependent_variable(psi);
      ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);

      std::cout << "Recorded data..." << std::endl;
      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
      std::cout << "v_ad: " << v_ad << std::endl;
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
      ad_helper.activate_recorded_tape(tape_no);
      v *= 2.2;
      ad_helper.set_independent_variable(v, v_dof);

      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
    }

  // Compute the function value, gradient and hessian for the new evaluation
  // point
  psi = ad_helper.compute_value();
  ad_helper.compute_gradient(Dpsi);
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      ad_helper.compute_hessian(D2psi);
    }

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
  const Tensor<1, dim, ScalarNumberType> dpsi_dv =
    ad_helper.extract_gradient_component(Dpsi, v_dof);

  // Verify the result
  using func = FunctionsTestVector<dim, ScalarNumberType>;
  static const ScalarNumberType tol =
    1e5 * std::numeric_limits<ScalarNumberType>::epsilon();
  std::cout << "psi:              " << psi << std::endl;
  std::cout << "func::psi(v):     " << func::psi(v) << std::endl;
  Assert(std::abs(psi - func::psi(v)) < tol,
         ExcMessage("No match for function value."));
  std::cout << "dpsi_dv:              " << dpsi_dv << std::endl;
  std::cout << "func::dpsi_dv(v):     " << func::dpsi_dv(v) << std::endl;
  Assert(std::abs((dpsi_dv - func::dpsi_dv(v)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      const Tensor<2, dim, ScalarNumberType> d2psi_dv_dv =
        ad_helper.extract_hessian_component(D2psi, v_dof, v_dof);
      std::cout << "d2psi_dv_dv:          " << d2psi_dv_dv << std::endl;
      std::cout << "func::d2psi_dv_dv(v): " << func::d2psi_dv_dv(v)
                << std::endl;
      Assert(std::abs((d2psi_dv_dv - func::d2psi_dv_dv(v)).norm()) < tol,
             ExcMessage("No match for second derivative."));
    }

  std::cout << std::endl << std::endl;
}
