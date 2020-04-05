// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Header file:
// Evaluation of a coupled system (vector + scalar components)
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
struct FunctionsTestVectorScalarCoupled
{
  static NumberType
  psi(const Tensor<1, dim, NumberType> &v, const NumberType &s)
  {
    return (v * v) * std::pow(s, 3);
  };

  static Tensor<1, dim, NumberType>
  dpsi_dv(const Tensor<1, dim, NumberType> &v, const NumberType &s)
  {
    return 2.0 * v * std::pow(s, 3);
  };

  static NumberType
  dpsi_ds(const Tensor<1, dim, NumberType> &v, const NumberType &s)
  {
    return 3.0 * (v * v) * std::pow(s, 2);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dv_dv(const Tensor<1, dim, NumberType> &v, const NumberType &s)
  {
    static const Tensor<2, dim, NumberType> I(
      unit_symmetric_tensor<dim, NumberType>());
    return 2.0 * I * std::pow(s, 3);
  };

  static Tensor<1, dim, NumberType>
  d2psi_ds_dv(const Tensor<1, dim, NumberType> &v, const NumberType &s)
  {
    return 6.0 * v * std::pow(s, 2);
  };

  static Tensor<1, dim, NumberType>
  d2psi_dv_ds(const Tensor<1, dim, NumberType> &v, const NumberType &s)
  {
    return d2psi_ds_dv(v, s);
  };

  static NumberType
  d2psi_ds_ds(const Tensor<1, dim, NumberType> &v, const NumberType &s)
  {
    return 6.0 * (v * v) * std::pow(s, 1);
  };
};

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_vector_scalar_coupled()
{
  typedef AD::ScalarFunction<dim, ad_type_code, number_t> ADHelper;
  typedef typename ADHelper::ad_type                      ADNumberType;
  typedef typename ADHelper::scalar_type                  ScalarNumberType;

  std::cout << "*** Test variables: Vector + Scalar (coupled), "
            << "dim = " << Utilities::to_string(dim) << ", "
            << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  // Values computed from the AD energy function
  ScalarNumberType             psi;
  Vector<ScalarNumberType>     Dpsi;
  FullMatrix<ScalarNumberType> D2psi;

  // Function and its derivatives
  typedef FunctionsTestVectorScalarCoupled<dim, ADNumberType> func_ad;

  // Setup the variable components and choose a value at which to
  // evaluate the tape
  const FEValuesExtractors::Vector v_dof(0);
  const FEValuesExtractors::Scalar s_dof(
    Tensor<1, dim>::n_independent_components);
  const unsigned int n_AD_components = dim + 1;
  ADHelper           ad_helper(n_AD_components);
  ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                     // default values

  ScalarNumberType                 s = 3.0;
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
      ad_helper.register_independent_variable(s, s_dof);

      const Tensor<1, dim, ADNumberType> v_ad =
        ad_helper.get_sensitive_variables(v_dof);
      const ADNumberType s_ad = ad_helper.get_sensitive_variables(s_dof);

      const ADNumberType psi(func_ad::psi(v_ad, s_ad));

      ad_helper.register_dependent_variable(psi);
      ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);

      std::cout << "Taped data..." << std::endl;
      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
      std::cout << "v_ad: " << v_ad << std::endl;
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
      ad_helper.activate_recorded_tape(tape_no);
      s = 2.75;
      v *= 1.3;
      ad_helper.set_independent_variable(v, v_dof);
      ad_helper.set_independent_variable(s, s_dof);
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
  const Tensor<1, dim, ScalarNumberType> dpsi_dv =
    ad_helper.extract_gradient_component(Dpsi, v_dof);
  const ScalarNumberType dpsi_ds =
    ad_helper.extract_gradient_component(Dpsi, s_dof);
  std::cout << "extracted Dpsi (v): " << dpsi_dv << "\n"
            << "extracted Dpsi (s): " << dpsi_ds << "\n";

  // Verify the result
  typedef FunctionsTestVectorScalarCoupled<dim, ScalarNumberType> func;
  static const ScalarNumberType                                   tol =
    1e5 * std::numeric_limits<ScalarNumberType>::epsilon();
  Assert(std::abs(psi - func::psi(v, s)) < tol,
         ExcMessage("No match for function value."));
  Assert(std::abs((dpsi_dv - func::dpsi_dv(v, s)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  Assert(std::abs(dpsi_ds - func::dpsi_ds(v, s)) < tol,
         ExcMessage("No match for first derivative."));
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      const Tensor<2, dim, ScalarNumberType> d2psi_dv_dv =
        ad_helper.extract_hessian_component(D2psi, v_dof, v_dof);
      const Tensor<1, dim, ScalarNumberType> d2psi_ds_dv =
        ad_helper.extract_hessian_component(D2psi, s_dof, v_dof);
      const Tensor<1, dim, ScalarNumberType> d2psi_dv_ds =
        ad_helper.extract_hessian_component(D2psi, v_dof, s_dof);
      const ScalarNumberType d2psi_ds_ds =
        ad_helper.extract_hessian_component(D2psi, s_dof, s_dof);
      std::cout << "extracted D2psi (v,v): " << d2psi_dv_dv << "\n"
                << "extracted D2psi (v,s): " << d2psi_ds_dv << "\n"
                << "extracted D2psi (s,v): " << d2psi_dv_ds << "\n"
                << "extracted D2psi (s,s): " << d2psi_ds_ds << "\n"
                << std::endl;
      Assert(std::abs((d2psi_dv_dv - func::d2psi_dv_dv(v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_ds_dv - func::d2psi_ds_dv(v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv_ds - func::d2psi_dv_ds(v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs(d2psi_ds_ds - func::d2psi_ds_ds(v, s)) < tol,
             ExcMessage("No match for second derivative."));
    }
}
