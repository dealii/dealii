// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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
// Test to check that tensor functions both compile and produce the right
// result when differentiated using the various auto-differentiable number
// types: Eigenvalues and eignvectors

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/differentiation/ad.h>

#include "../tests.h"
#ifdef DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING
#  include <deal.II/base/symmetric_tensor.templates.h>
#endif

#include <iostream>

namespace AD = dealii::Differentiation::AD;

double
mu()
{
  return 100.0;
}


template <int dim, typename NumberType>
NumberType
J(const SymmetricTensor<2, dim, NumberType> &C)
{
  return std::sqrt(determinant(C));
}

template <int dim, typename NumberType>
struct IncompressibleNeoHookean
{
  // Incompressible Neo-Hookean material
  static NumberType
  psi(const SymmetricTensor<2, dim, NumberType> &C)
  {
    return 0.5 * mu() * (trace(C) - dim) - mu() * std::log(J(C));
  }

  static SymmetricTensor<2, dim, NumberType>
  dpsi_dC(const SymmetricTensor<2, dim, NumberType> &C)
  {
    const SymmetricTensor<2, dim, NumberType> I =
      unit_symmetric_tensor<dim, NumberType>();
    return 0.5 * mu() * (I - invert(C));
  }

  static SymmetricTensor<4, dim, NumberType>
  d2psi_dC_dC(const SymmetricTensor<2, dim, NumberType> &C)
  {
    const SymmetricTensor<2, dim, NumberType> C_inv = invert(C);

    SymmetricTensor<4, dim, NumberType> dC_inv_dC;
    for (unsigned int A = 0; A < dim; ++A)
      for (unsigned int B = A; B < dim; ++B)
        for (unsigned int C = 0; C < dim; ++C)
          for (unsigned int D = C; D < dim; ++D)
            dC_inv_dC[A][B][C][D] -=
              0.5 * (C_inv[A][C] * C_inv[B][D] + C_inv[A][D] * C_inv[B][C]);

    return -0.5 * mu() * dC_inv_dC;
  }
};

template <int dim, typename NumberType>
struct IncompressibleNeoHookeanPrincipalStretches
{
  //  // Incompressible Neo-Hookean material
  //  static NumberType
  //  psi (const std::array<std::pair<NumberType, Tensor<1,dim,NumberType>
  //  >,dim> eig_C)
  //  {
  //    NumberType psi = 0.0;
  //    NumberType J = 1.0;
  //    for (unsigned int d=0; d<dim; ++d)
  //      {
  //        const NumberType &lambda_squared = eig_C[d].first;
  //        psi += 0.5*mu()*(lambda_squared - 1.0);
  //        J *= std::sqrt(lambda_squared);
  //      }
  //    psi -= mu()*std::log(J);
  //    return psi;
  //  }

  static SymmetricTensor<2, dim, NumberType>
  dpsi_dC(const std::array<std::pair<NumberType, Tensor<1, dim, NumberType>>,
                           dim> eig_C)
  {
    SymmetricTensor<2, dim, NumberType> C_inv;
    for (unsigned int d = 0; d < dim; ++d)
      {
        const NumberType &                lambda_squared = eig_C[d].first;
        const Tensor<1, dim, NumberType> &N              = eig_C[d].second;
        C_inv += (1.0 / lambda_squared) * symmetrize(outer_product(N, N));
      }

    const SymmetricTensor<2, dim, NumberType> I =
      unit_symmetric_tensor<dim, NumberType>();
    return 0.5 * mu() * (I - C_inv);
  }
};

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_NH(const bool nontrivial_initial_values)
{
  typedef AD::ScalarFunction<dim, ad_type_code, number_t> ADHelper;
  typedef typename ADHelper::ad_type                      ADNumberType;
  typedef typename ADHelper::scalar_type                  ScalarNumberType;

  std::cout << "*** Standard definition of incompressible NeoHookean material, "
            << "dim = " << Utilities::to_string(dim) << ", "
            << "Type code: " << static_cast<int>(ad_type_code) << ", "
            << "Nontrivial initial values: " << std::boolalpha
            << nontrivial_initial_values << std::endl;

  // Values computed from the AD energy function
  ScalarNumberType             psi;
  Vector<ScalarNumberType>     Dpsi;
  FullMatrix<ScalarNumberType> D2psi;

  // Function and its derivatives
  typedef IncompressibleNeoHookean<dim, ADNumberType> func_ad;

  // Setup the variable components and choose a value at which to
  // evaluate the tape
  const FEValuesExtractors::SymmetricTensor<2> C_dof(0);
  const unsigned int                           n_AD_components =
    SymmetricTensor<2, dim>::n_independent_components;
  ADHelper ad_helper(n_AD_components);
  ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                     // default values

  SymmetricTensor<2, dim, ScalarNumberType> C =
    unit_symmetric_tensor<dim, ScalarNumberType>();
  if (nontrivial_initial_values)
    for (unsigned int i = 0; i < C.n_independent_components; ++i)
      C[C.unrolled_to_component_indices(i)] += 0.12 * (i + 0.02);

  const int  tape_no = 1;
  const bool is_recording =
    ad_helper.start_recording_operations(tape_no /*material_id*/,
                                         true /*overwrite_tape*/,
                                         true /*keep*/);
  if (is_recording == true)
    {
      ad_helper.register_independent_variable(C, C_dof);

      const SymmetricTensor<2, dim, ADNumberType> C_ad =
        ad_helper.get_sensitive_variables(C_dof);

      const ADNumberType psi(func_ad::psi(C_ad));

      ad_helper.register_dependent_variable(psi);
      ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);

      std::cout << "Recorded data..." << std::endl;
      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
      std::cout << "C_ad: " << C_ad << std::endl;
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
      C *= 1.15;
      ad_helper.set_independent_variable(C, C_dof);

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
  const SymmetricTensor<2, dim, ScalarNumberType> dpsi_dC =
    ad_helper.extract_gradient_component(Dpsi, C_dof);

  // Verify the result
  typedef IncompressibleNeoHookean<dim, ScalarNumberType> func;
  static const ScalarNumberType                           tol =
    1e5 * std::numeric_limits<ScalarNumberType>::epsilon();
  std::cout << "psi:              " << psi << std::endl;
  std::cout << "func::psi(C):     " << func::psi(C) << std::endl;
  Assert(std::abs(psi - func::psi(C)) < tol,
         ExcMessage("No match for function value."));
  std::cout << "dpsi_dC:              " << dpsi_dC << std::endl;
  std::cout << "func::dpsi_dC(C):     " << func::dpsi_dC(C) << std::endl;
  Assert(std::abs((dpsi_dC - func::dpsi_dC(C)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      const Tensor<4, dim, ScalarNumberType> d2psi_dC_dC =
        ad_helper.extract_hessian_component(D2psi, C_dof, C_dof);
      std::cout << "d2psi_dC_dC:          " << d2psi_dC_dC << std::endl;
      std::cout << "func::d2psi_dC_dC(C): " << func::d2psi_dC_dC(C)
                << std::endl;
      Assert(std::abs((d2psi_dC_dC - func::d2psi_dC_dC(C)).norm()) < tol,
             ExcMessage("No match for second derivative."));
    }

  std::cout << std::endl << std::endl;
}



template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_NH_eigen_stress(const enum SymmetricTensorEigenvectorMethod method,
                     const bool nontrivial_initial_values)
{
  typedef AD::VectorFunction<dim, ad_type_code, number_t> ADHelper;
  typedef typename ADHelper::ad_type                      ADNumberType;
  typedef typename ADHelper::scalar_type                  ScalarNumberType;

  std::cout
    << "*** Principal stretch definition of incompressible NeoHookean material (from stress), "
    << "dim = " << Utilities::to_string(dim) << ", "
    << "Type code: " << static_cast<int>(ad_type_code) << ", "
    << "Eig method: " << static_cast<int>(method) << ", "
    << "Nontrivial initial values: " << std::boolalpha
    << nontrivial_initial_values << std::endl;

  // Values computed from the AD energy function
  // ScalarNumberType psi;
  Vector<ScalarNumberType>     Dpsi;
  FullMatrix<ScalarNumberType> D2psi;

  // Function and its derivatives
  typedef IncompressibleNeoHookeanPrincipalStretches<dim, ADNumberType> func_ad;

  // Setup the variable components and choose a value at which to
  // evaluate the tape
  const FEValuesExtractors::SymmetricTensor<2> C_dof(0);
  const unsigned int                           n_AD_independent_components =
    SymmetricTensor<2, dim>::n_independent_components;
  const unsigned int n_AD_dependent_components =
    SymmetricTensor<2, dim>::n_independent_components;
  ADHelper ad_helper(n_AD_independent_components, n_AD_dependent_components);
  ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                     // default values

  SymmetricTensor<2, dim, ScalarNumberType> C =
    unit_symmetric_tensor<dim, ScalarNumberType>();
  if (nontrivial_initial_values)
    for (unsigned int i = 0; i < C.n_independent_components; ++i)
      C[C.unrolled_to_component_indices(i)] += 0.12 * (i + 0.02);

  const int  tape_no = 1;
  const bool is_recording =
    ad_helper.start_recording_operations(tape_no /*material_id*/,
                                         true /*overwrite_tape*/,
                                         true /*keep*/);
  if (is_recording == true)
    {
      ad_helper.register_independent_variable(C, C_dof);

      const SymmetricTensor<2, dim, ADNumberType> C_ad =
        ad_helper.get_sensitive_variables(C_dof);
      const auto eig_C_ad = eigenvectors(C_ad, method);

      for (unsigned int d = 0; d < dim; ++d)
        std::cout << "  Direction: " << d
                  << "  Eigenvalue: " << eig_C_ad[d].first
                  << "  Eigenvector: " << eig_C_ad[d].second << std::endl;

      const SymmetricTensor<2, dim, ADNumberType> dpsi_dC(
        func_ad::dpsi_dC(eig_C_ad));

      ad_helper.register_dependent_variable(dpsi_dC, C_dof);
      ad_helper.stop_recording_operations(false /*write_tapes_to_file*/);

      std::cout << "Recorded data..." << std::endl;
      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
      std::cout << "C_ad: " << C_ad << std::endl;
      std::cout << "dpsi_dC: " << dpsi_dC << std::endl;
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
      C *= 1.15;
      ad_helper.set_independent_variable(C, C_dof);

      std::cout << "independent variable values: " << std::flush;
      ad_helper.print_values(std::cout);
    }

  // Compute the function value, gradient and hessian for the new evaluation
  // point
  // psi = ad_helper.compute_value();
  ad_helper.compute_values(Dpsi);
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1)
    {
      ad_helper.compute_jacobian(D2psi);
    }

  // Output the full stored function, gradient vector and hessian matrix
  // std::cout << "psi: " << psi << std::endl;
  std::cout << "Dpsi: \n";
  Dpsi.print(std::cout);
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1)
    {
      std::cout << "D2psi: \n";
      D2psi.print_formatted(std::cout, 3, true, 0, "0.0");
    }

  // Extract components of the solution
  const SymmetricTensor<2, dim, ScalarNumberType> dpsi_dC =
    ad_helper.extract_value_component(Dpsi, C_dof);

  // Verify the result
  typedef IncompressibleNeoHookean<dim, ScalarNumberType> func;
  static const ScalarNumberType                           tol_val =
    (nontrivial_initial_values ?
       1e-4 :
       1e6 * std::numeric_limits<ScalarNumberType>::epsilon());
  static const ScalarNumberType tol_jac =
    (nontrivial_initial_values ?
       2.5e-3 :
       1e6 * std::numeric_limits<ScalarNumberType>::epsilon());
  //  std::cout << "psi:              " << psi << std::endl;
  //  std::cout << "func::psi(C):     " << func::psi(C) << std::endl;
  //  Assert(std::abs(psi - func::psi(C)) < tol, ExcMessage("No match for
  //  function value."));
  std::cout << "dpsi_dC:              " << dpsi_dC << std::endl;
  std::cout << "func::dpsi_dC(C):     " << func::dpsi_dC(C) << std::endl;
  //  std::cout << "DIFF NORM: " << std::abs((dpsi_dC -
  //  func::dpsi_dC(C)).norm()) << std::endl;
  Assert(std::abs((dpsi_dC - func::dpsi_dC(C)).norm()) < tol_val,
         ExcMessage("No match for first derivative."));
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 1)
    {
      const SymmetricTensor<4, dim, ScalarNumberType> d2psi_dC_dC =
        ad_helper.extract_jacobian_component(D2psi, C_dof, C_dof);
      std::cout << "d2psi_dC_dC:          " << d2psi_dC_dC << std::endl;
      std::cout << "func::d2psi_dC_dC(C): " << func::d2psi_dC_dC(C)
                << std::endl;
      //      std::cout << "DIFF: " << (d2psi_dC_dC - func::d2psi_dC_dC(C)) <<
      //      std::endl; std::cout << "DIFF NORM: " << std::abs((d2psi_dC_dC -
      //      func::d2psi_dC_dC(C)).norm()) << std::endl;
      Assert(std::abs((d2psi_dC_dC - func::d2psi_dC_dC(C)).norm()) < tol_jac,
             ExcMessage("No match for second derivative."));
    }

  std::cout << std::endl << std::endl;
}



template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_symmetric_tensor()
{
  // Then test the functions that compute the eigenvalues/vectors:
  // Note: Derived from a stress, both the eigenvalues and eigenvectors are used

  // Non-trivial initial values (unequal eigenvalues):
  {
    const bool nontrivial_initial_values = true;

    // First verify that all manual calculations are correct
    test_NH<dim, number_t, ad_type_code>(nontrivial_initial_values);

    test_NH_eigen_stress<dim, number_t, ad_type_code>(
      SymmetricTensorEigenvectorMethod::hybrid, nontrivial_initial_values);
    test_NH_eigen_stress<dim, number_t, ad_type_code>(
      SymmetricTensorEigenvectorMethod::ql_implicit_shifts,
      nontrivial_initial_values);
    test_NH_eigen_stress<dim, number_t, ad_type_code>(
      SymmetricTensorEigenvectorMethod::jacobi, nontrivial_initial_values);
  }
  // Trivial initial values (equal eigenvalues):
  {
    const bool nontrivial_initial_values = false;

    // First verify that all manual calculations are correct
    test_NH<dim, number_t, ad_type_code>(nontrivial_initial_values);

    //    test_NH_eigen_stress<dim,number_t,ad_type_code>(SymmetricTensorEigenvectorMethod::hybrid,nontrivial_initial_values);
    //    // This will never work.
    test_NH_eigen_stress<dim, number_t, ad_type_code>(
      SymmetricTensorEigenvectorMethod::ql_implicit_shifts,
      nontrivial_initial_values);
    test_NH_eigen_stress<dim, number_t, ad_type_code>(
      SymmetricTensorEigenvectorMethod::jacobi, nontrivial_initial_values);
  }
}
