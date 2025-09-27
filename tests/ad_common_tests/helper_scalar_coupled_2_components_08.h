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
// Evaluation of a coupled system (symmetric tensor + symmetric tensor
// components) using a helper class

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
struct FunctionsTestSymmetricTensorSymmetricTensorCoupled
{
  static NumberType
  det_t(const SymmetricTensor<2, dim, NumberType> &t)
  {
    return determinant(t);
  }

  static SymmetricTensor<2, dim, NumberType>
  ddet_t_dt(const SymmetricTensor<2, dim, NumberType> &t)
  {
    const SymmetricTensor<2, dim, NumberType> t_inv =
      symmetrize(invert(static_cast<Tensor<2, dim, NumberType>>(t)));
    return det_t(t) * t_inv;
  }

  static SymmetricTensor<4, dim, NumberType>
  d2det_t_dt_dt(const SymmetricTensor<2, dim, NumberType> &t)
  {
    const SymmetricTensor<2, dim, NumberType> t_inv =
      symmetrize(invert(static_cast<Tensor<2, dim, NumberType>>(t)));
    SymmetricTensor<4, dim, NumberType> dt_inv_dt;
    // https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Tensors#Derivative_of_the_determinant_of_a_tensor
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = k; l < dim; ++l)
            dt_inv_dt[i][j][k][l] =
              -0.5 * (t_inv[i][k] * t_inv[j][l] + t_inv[i][l] * t_inv[j][k]);

    return det_t(t) * outer_product(t_inv, t_inv) + det_t(t) * dt_inv_dt;
  }

  // -------

  static NumberType
  psi(const SymmetricTensor<2, dim, NumberType> &t1,
      const SymmetricTensor<2, dim, NumberType> &t2)
  {
    return pow(det_t(t1), 2) * pow(det_t(t2), 3);
  };

  static SymmetricTensor<2, dim, NumberType>
  dpsi_dt1(const SymmetricTensor<2, dim, NumberType> &t1,
           const SymmetricTensor<2, dim, NumberType> &t2)
  {
    return SymmetricTensor<2, dim, NumberType>(
      2.0 * pow(det_t(t1), 1) * ddet_t_dt(t1) * pow(det_t(t2), 3));
  };

  static SymmetricTensor<2, dim, NumberType>
  dpsi_dt2(const SymmetricTensor<2, dim, NumberType> &t1,
           const SymmetricTensor<2, dim, NumberType> &t2)
  {
    return SymmetricTensor<2, dim, NumberType>(
      pow(det_t(t1), 2) * 3.0 * pow(det_t(t2), 2) * ddet_t_dt(t2));
  };

  static SymmetricTensor<4, dim, NumberType>
  d2psi_dt1_dt1(const SymmetricTensor<2, dim, NumberType> &t1,
                const SymmetricTensor<2, dim, NumberType> &t2)
  {
    return SymmetricTensor<4, dim, NumberType>(
      2.0 * pow(det_t(t2), 3) *
      (pow(det_t(t1), 0) * outer_product(ddet_t_dt(t1), ddet_t_dt(t1)) +
       pow(det_t(t1), 1) * d2det_t_dt_dt(t1)));
  };

  static SymmetricTensor<4, dim, NumberType>
  d2psi_dt2_dt1(const SymmetricTensor<2, dim, NumberType> &t1,
                const SymmetricTensor<2, dim, NumberType> &t2)
  {
    return SymmetricTensor<4, dim, NumberType>(
      2.0 * pow(det_t(t1), 1) * 3.0 * pow(det_t(t2), 2) *
      outer_product(ddet_t_dt(t1), ddet_t_dt(t2)));
  };

  static SymmetricTensor<4, dim, NumberType>
  d2psi_dt1_dt2(const SymmetricTensor<2, dim, NumberType> &t1,
                const SymmetricTensor<2, dim, NumberType> &t2)
  {
    return SymmetricTensor<4, dim, NumberType>(
      2.0 * pow(det_t(t1), 1) * 3.0 * pow(det_t(t2), 2) *
      outer_product(ddet_t_dt(t2), ddet_t_dt(t1)));
  };

  static SymmetricTensor<4, dim, NumberType>
  d2psi_dt2_dt2(const SymmetricTensor<2, dim, NumberType> &t1,
                const SymmetricTensor<2, dim, NumberType> &t2)
  {
    return SymmetricTensor<4, dim, NumberType>(
      pow(det_t(t1), 2) * 3.0 *
      (2.0 * pow(det_t(t2), 1) * outer_product(ddet_t_dt(t2), ddet_t_dt(t2)) +
       pow(det_t(t2), 2) * d2det_t_dt_dt(t2)));
  };
};

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_symmetric_tensor_symmetric_tensor_coupled()
{
  using ADHelper         = AD::ScalarFunction<dim, ad_type_code, number_t>;
  using ADNumberType     = typename ADHelper::ad_type;
  using ScalarNumberType = typename ADHelper::scalar_type;

  std::cout
    << "*** Test variables: SymmetricTensor + SymmetricTensor (coupled), "
    << "dim = " << Utilities::to_string(dim) << ", "
    << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  // Values computed from the AD energy function
  ScalarNumberType             psi;
  Vector<ScalarNumberType>     Dpsi;
  FullMatrix<ScalarNumberType> D2psi;

  // Function and its derivatives
  using func_ad =
    FunctionsTestSymmetricTensorSymmetricTensorCoupled<dim, ADNumberType>;

  // Setup the variable components and choose a value at which to
  // evaluate the tape
  const FEValuesExtractors::SymmetricTensor<2> t1_dof(0);
  const FEValuesExtractors::SymmetricTensor<2> t2_dof(
    SymmetricTensor<2, dim>::n_independent_components);
  const unsigned int n_AD_components =
    SymmetricTensor<2, dim>::n_independent_components +
    SymmetricTensor<2, dim>::n_independent_components;
  ADHelper ad_helper(n_AD_components);
  ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                     // default values

  SymmetricTensor<2, dim, ScalarNumberType> t1 =
    unit_symmetric_tensor<dim, ScalarNumberType>();
  for (unsigned int i = 0; i < t1.n_independent_components; ++i)
    t1[t1.unrolled_to_component_indices(i)] += 0.135 * (i + 0.195);
  SymmetricTensor<2, dim, ScalarNumberType> t2 =
    unit_symmetric_tensor<dim, ScalarNumberType>();
  for (unsigned int i = 0; i < t2.n_independent_components; ++i)
    t2[t2.unrolled_to_component_indices(i)] += 0.18 * (i + 0.09);

  const int  tape_no = 1;
  const bool is_recording =
    ad_helper.start_recording_operations(tape_no /*material_id*/,
                                         true /*overwrite_tape*/,
                                         true /*keep*/);
  if (is_recording == true)
    {
      ad_helper.register_independent_variable(t1, t1_dof);
      ad_helper.register_independent_variable(t2, t2_dof);

      const SymmetricTensor<2, dim, ADNumberType> t1_ad =
        ad_helper.get_sensitive_variables(t1_dof);
      const SymmetricTensor<2, dim, ADNumberType> t2_ad =
        ad_helper.get_sensitive_variables(t2_dof);

      const ADNumberType psi(func_ad::psi(t1_ad, t2_ad));

      ad_helper.register_dependent_variable(psi);
      ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);

      std::cout << "Recorded data..." << std::endl;
      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
      std::cout << "t1_ad: " << t1_ad << std::endl;
      std::cout << "t2_ad: " << t2_ad << std::endl;
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
      t1 *= 1.17;
      t2 *= 0.92;
      ad_helper.set_independent_variable(t1, t1_dof);
      ad_helper.set_independent_variable(t2, t2_dof);
    }

  std::cout << "independent variable values: " << std::flush;
  ad_helper.print_values(std::cout);

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
  const SymmetricTensor<2, dim, ScalarNumberType> dpsi_dt1 =
    ad_helper.extract_gradient_component(Dpsi, t1_dof);
  const SymmetricTensor<2, dim, ScalarNumberType> dpsi_dt2 =
    ad_helper.extract_gradient_component(Dpsi, t2_dof);
  std::cout << "extracted Dpsi (t1): " << dpsi_dt1 << "\n"
            << "extracted Dpsi (t2): " << dpsi_dt2 << "\n";

  // Verify the result
  using func =
    FunctionsTestSymmetricTensorSymmetricTensorCoupled<dim, ScalarNumberType>;
  static const ScalarNumberType tol =
    1e5 * std::numeric_limits<ScalarNumberType>::epsilon();
  std::cout << "dpsi_dt1:            " << dpsi_dt1 << std::endl;
  std::cout << "func::dpsi_dt1(t1,t2): " << func::dpsi_dt1(t1, t2) << std::endl;
  std::cout << "diff: " << std::abs((dpsi_dt1 - func::dpsi_dt1(t1, t2)).norm())
            << std::endl;
  std::cout << "dpsi_dt2:            " << dpsi_dt2 << std::endl;
  std::cout << "func::dpsi_dt2(t1,t2): " << func::dpsi_dt2(t1, t2) << std::endl;
  std::cout << "diff: " << std::abs((dpsi_dt2 - func::dpsi_dt2(t1, t2)).norm())
            << std::endl;
  Assert(std::abs(psi - func::psi(t1, t2)) < tol,
         ExcMessage("No match for function value."));
  Assert(std::abs((dpsi_dt1 - func::dpsi_dt1(t1, t2)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  Assert(std::abs((dpsi_dt2 - func::dpsi_dt2(t1, t2)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      const SymmetricTensor<4, dim, ScalarNumberType> d2psi_dt1_dt1 =
        ad_helper.extract_hessian_component(D2psi, t1_dof, t1_dof);
      const SymmetricTensor<4, dim, ScalarNumberType> d2psi_dt2_dt1 =
        ad_helper.extract_hessian_component(D2psi, t1_dof, t2_dof);
      const SymmetricTensor<4, dim, ScalarNumberType> d2psi_dt1_dt2 =
        ad_helper.extract_hessian_component(D2psi, t2_dof, t1_dof);
      const SymmetricTensor<4, dim, ScalarNumberType> d2psi_dt2_dt2 =
        ad_helper.extract_hessian_component(D2psi, t2_dof, t2_dof);
      std::cout << "extracted Dpsi (t1): " << dpsi_dt1 << "\n"
                << "extracted Dpsi (t2): " << dpsi_dt2 << "\n"
                << "extracted D2psi (t1,t1): " << d2psi_dt1_dt1 << "\n"
                << "extracted D2psi (t1,t2): " << d2psi_dt2_dt1 << "\n"
                << "extracted D2psi (t2,t1): " << d2psi_dt1_dt2 << "\n"
                << "extracted D2psi (t2,t2): " << d2psi_dt2_dt2 << "\n"
                << std::endl;
      std::cout << "d2psi_dt1_dt1:            " << d2psi_dt1_dt1 << std::endl;
      std::cout << "func::d2psi_dt1_dt1(t1,t2): " << func::d2psi_dt1_dt1(t1, t2)
                << std::endl;
      std::cout << "diff: "
                << std::abs(
                     (d2psi_dt1_dt1 - func::d2psi_dt1_dt1(t1, t2)).norm())
                << std::endl;
      std::cout << "d2psi_dt2_dt1:            " << d2psi_dt2_dt1 << std::endl;
      std::cout << "func::d2psi_dt2_dt1(t1,t2): " << func::d2psi_dt2_dt1(t1, t2)
                << std::endl;
      std::cout << "diff: "
                << std::abs(
                     (d2psi_dt2_dt1 - func::d2psi_dt2_dt1(t1, t2)).norm())
                << std::endl;
      std::cout << "d2psi_dt1_dt2:            " << d2psi_dt1_dt2 << std::endl;
      std::cout << "func::d2psi_dt1_dt2(t1,t2): " << func::d2psi_dt1_dt2(t1, t2)
                << std::endl;
      std::cout << "diff: "
                << std::abs(
                     (d2psi_dt1_dt2 - func::d2psi_dt1_dt2(t1, t2)).norm())
                << std::endl;
      std::cout << "d2psi_dt2_dt2:            " << d2psi_dt2_dt2 << std::endl;
      std::cout << "func::d2psi_dt2_dt2(t1,t2): " << func::d2psi_dt2_dt2(t1, t2)
                << std::endl;
      std::cout << "diff: "
                << std::abs(
                     (d2psi_dt2_dt2 - func::d2psi_dt2_dt2(t1, t2)).norm())
                << std::endl;
      Assert(std::abs((d2psi_dt1_dt1 - func::d2psi_dt1_dt1(t1, t2)).norm()) <
               tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dt2_dt1 - func::d2psi_dt2_dt1(t1, t2)).norm()) <
               tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dt1_dt2 - func::d2psi_dt1_dt2(t1, t2)).norm()) <
               tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dt2_dt2 - func::d2psi_dt2_dt2(t1, t2)).norm()) <
               tol,
             ExcMessage("No match for second derivative."));
    }
}
