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
// Evaluation of a coupled system (vector + vector components)
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
struct FunctionsTestVectorVectorCoupled
{
  static NumberType
  psi(const Tensor<1, dim, NumberType> &v1,
      const Tensor<1, dim, NumberType> &v2)
  {
    return std::pow(v1 * v1, 2) * std::pow(v2 * v2, 3);
  };

  static Tensor<1, dim, NumberType>
  dpsi_dv1(const Tensor<1, dim, NumberType> &v1,
           const Tensor<1, dim, NumberType> &v2)
  {
    return 2.0 * std::pow(v1 * v1, 1) * 2.0 * v1 * std::pow(v2 * v2, 3);
  };

  static Tensor<1, dim, NumberType>
  dpsi_dv2(const Tensor<1, dim, NumberType> &v1,
           const Tensor<1, dim, NumberType> &v2)
  {
    return std::pow(v1 * v1, 2) * 3.0 * std::pow(v2 * v2, 2) * 2.0 * v2;
  };

  static Tensor<2, dim, NumberType>
  d2psi_dv1_dv1(const Tensor<1, dim, NumberType> &v1,
                const Tensor<1, dim, NumberType> &v2)
  {
    const Tensor<2, dim, NumberType> I(
      unit_symmetric_tensor<dim, NumberType>());
    return 2.0 * 2.0 * std::pow(v2 * v2, 3) *
           (pow(v1 * v1, 0) * 2.0 * outer_product(v1, v1) +
            std::pow(v1 * v1, 1) * I);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dv2_dv1(const Tensor<1, dim, NumberType> &v1,
                const Tensor<1, dim, NumberType> &v2)
  {
    // Note: This is not symmetric, and this is why we
    // don't set the return type for hessian extractor (vector,vector)
    // operations as SymmetricTensor.
    return (2.0 * std::pow(v1 * v1, 1) * 2.0) *
           (3.0 * std::pow(v2 * v2, 2) * 2.0) * outer_product(v1, v2);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dv1_dv2(const Tensor<1, dim, NumberType> &v1,
                const Tensor<1, dim, NumberType> &v2)
  {
    return (2.0 * std::pow(v1 * v1, 1) * 2.0) *
           (3.0 * std::pow(v2 * v2, 2) * 2.0) * outer_product(v2, v1);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dv2_dv2(const Tensor<1, dim, NumberType> &v1,
                const Tensor<1, dim, NumberType> &v2)
  {
    const Tensor<2, dim, NumberType> I(
      unit_symmetric_tensor<dim, NumberType>());
    // return std::pow(v1*v1,2)*3.0*
    //        ( 2.0*std::pow(v2*v2,1)*2.0*outer_product(v2,v2) +
    //        std::pow(v2*v2,2)*2.0*I);
    return std::pow(v1 * v1, 2) * 3.0 * 2.0 *
           (2.0 * std::pow(v2 * v2, 1) * 2.0 * outer_product(v2, v2) +
            std::pow(v2 * v2, 2) * I);
  };
};

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_vector_vector_coupled()
{
  typedef AD::ScalarFunction<dim, ad_type_code, number_t> ADHelper;
  typedef typename ADHelper::ad_type                      ADNumberType;
  typedef typename ADHelper::scalar_type                  ScalarNumberType;

  std::cout << "*** Test variables: Vector + Vector (coupled), "
            << "dim = " << Utilities::to_string(dim) << ", "
            << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  // Values computed from the AD energy function
  ScalarNumberType             psi;
  Vector<ScalarNumberType>     Dpsi;
  FullMatrix<ScalarNumberType> D2psi;

  // Function and its derivatives
  typedef FunctionsTestVectorVectorCoupled<dim, ADNumberType> func_ad;

  const FEValuesExtractors::Vector v1_dof(0);
  const FEValuesExtractors::Vector v2_dof(
    Tensor<1, dim>::n_independent_components);
  const unsigned int n_AD_components = dim + dim;
  ADHelper           ad_helper(n_AD_components);
  ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                     // default values

  Tensor<1, dim, ScalarNumberType> v1, v2;
  for (unsigned int i = 0; i < dim; ++i)
    {
      v1[i] = 0.2 * (1.0 + i);
      v2[i] = 0.1 * (2.0 + i * 1.7);
    }

  const int  tape_no = 1;
  const bool is_recording =
    ad_helper.start_recording_operations(tape_no /*material_id*/,
                                         true /*overwrite_tape*/,
                                         true /*keep*/);
  if (is_recording == true)
    {
      ad_helper.register_independent_variable(v1, v1_dof);
      ad_helper.register_independent_variable(v2, v2_dof);

      const Tensor<1, dim, ADNumberType> v1_ad =
        ad_helper.get_sensitive_variables(v1_dof);
      const Tensor<1, dim, ADNumberType> v2_ad =
        ad_helper.get_sensitive_variables(v2_dof);

      const ADNumberType psi(func_ad::psi(v1_ad, v2_ad));

      ad_helper.register_dependent_variable(psi);
      ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);

      std::cout << "Recorded data..." << std::endl;
      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
      std::cout << "v1_ad: " << v1_ad << std::endl;
      std::cout << "v2_ad: " << v2_ad << std::endl;
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
      v1 *= 1.7;
      v2 *= 0.25;
      ad_helper.set_independent_variable(v1, v1_dof);
      ad_helper.set_independent_variable(v2, v2_dof);
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
  const Tensor<1, dim, ScalarNumberType> dpsi_dv1 =
    ad_helper.extract_gradient_component(Dpsi, v1_dof);
  const Tensor<1, dim, ScalarNumberType> dpsi_dv2 =
    ad_helper.extract_gradient_component(Dpsi, v2_dof);
  std::cout << "extracted Dpsi (v1): " << dpsi_dv1 << "\n"
            << "extracted Dpsi (v2): " << dpsi_dv2 << "\n";

  // Verify the result
  typedef FunctionsTestVectorVectorCoupled<dim, ScalarNumberType> func;
  static const ScalarNumberType                                   tol =
    1e5 * std::numeric_limits<ScalarNumberType>::epsilon();

  Assert(std::abs(psi - func::psi(v1, v2)) < tol,
         ExcMessage("No match for function value."));
  Assert(std::abs((dpsi_dv1 - func::dpsi_dv1(v1, v2)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  Assert(std::abs((dpsi_dv2 - func::dpsi_dv2(v1, v2)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      const Tensor<2, dim, ScalarNumberType> d2psi_dv1_dv1 =
        ad_helper.extract_hessian_component(D2psi, v1_dof, v1_dof);
      const Tensor<2, dim, ScalarNumberType> d2psi_dv2_dv1 =
        ad_helper.extract_hessian_component(D2psi, v1_dof, v2_dof);
      const Tensor<2, dim, ScalarNumberType> d2psi_dv1_dv2 =
        ad_helper.extract_hessian_component(D2psi, v2_dof, v1_dof);
      const Tensor<2, dim, ScalarNumberType> d2psi_dv2_dv2 =
        ad_helper.extract_hessian_component(D2psi, v2_dof, v2_dof);
      std::cout << "extracted D2psi (v1,v1): " << d2psi_dv1_dv1 << "\n"
                << "extracted D2psi (v1,v2): " << d2psi_dv2_dv1 << "\n"
                << "extracted D2psi (v2,v1): " << d2psi_dv1_dv2 << "\n"
                << "extracted D2psi (v2,v2): " << d2psi_dv2_dv2 << "\n"
                << std::endl;
      Assert(std::abs((d2psi_dv1_dv1 - func::d2psi_dv1_dv1(v1, v2)).norm()) <
               tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv2_dv1 - func::d2psi_dv2_dv1(v1, v2)).norm()) <
               tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv1_dv2 - func::d2psi_dv1_dv2(v1, v2)).norm()) <
               tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv2_dv2 - func::d2psi_dv2_dv2(v1, v2)).norm()) <
               tol,
             ExcMessage("No match for second derivative."));
    }
}
