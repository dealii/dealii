// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// Test invariants: Coupled isotropic

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/physics/invariants.h>

#include <functional>

#include "../tests.h"


constexpr Differentiation::AD::NumberTypes ADTypeCode =
  Differentiation::AD::NumberTypes::sacado_dfad_dfad;
template <int dim>
using ADHelper = Differentiation::AD::ScalarFunction<dim, ADTypeCode, double>;
using ADNumberType = typename ADHelper<1>::ad_type;
template <int dim>
using psi_function_type =
  std::function<ADNumberType(const SymmetricTensor<2, dim, ADNumberType> &,
                             const Tensor<1, dim, ADNumberType> &)>;

template <int dim>
struct Values
{
  double                  Psi;
  SymmetricTensor<2, dim> dPsi_dC;
  Tensor<1, dim>          dPsi_dH;
  SymmetricTensor<4, dim> d2Psi_dC_dC;
  Tensor<3, dim>          d2Psi_dC_dH;
  Tensor<3, dim>          d2Psi_dH_dC;
  SymmetricTensor<2, dim> d2Psi_dH_dH;
};

template <int dim>
Values<dim>
compute_derivatives_using_AD(const SymmetricTensor<2, dim> &C,
                             const Tensor<1, dim> &         H,
                             const psi_function_type<dim> & get_psi_ad)
{
  const FEValuesExtractors::SymmetricTensor<2> C_components(0);
  const FEValuesExtractors::Vector H_components(C.n_independent_components);
  const unsigned int               n_independent_variables =
    C.n_independent_components + H.n_independent_components;
  ADHelper<dim> ad_helper(n_independent_variables);

  ad_helper.register_independent_variable(C, C_components);
  ad_helper.register_independent_variable(H, H_components);
  const SymmetricTensor<2, dim, ADNumberType> C_ad =
    ad_helper.get_sensitive_variables(C_components);
  const Tensor<1, dim, ADNumberType> H_ad =
    ad_helper.get_sensitive_variables(H_components);

  const ADNumberType psi_ad = get_psi_ad(C_ad, H_ad);
  ad_helper.register_dependent_variable(psi_ad);

  Vector<double>     Dpsi(ad_helper.n_independent_variables());
  FullMatrix<double> D2psi(ad_helper.n_independent_variables(),
                           ad_helper.n_independent_variables());

  ad_helper.compute_gradient(Dpsi);
  ad_helper.compute_hessian(D2psi);

  return {
    ad_helper.compute_value(),
    ad_helper.extract_gradient_component(Dpsi, C_components),
    ad_helper.extract_gradient_component(Dpsi, H_components),
    ad_helper.extract_hessian_component(D2psi, C_components, C_components),
    ad_helper.extract_hessian_component(D2psi, C_components, H_components),
    ad_helper.extract_hessian_component(D2psi, H_components, C_components),
    symmetrize(
      ad_helper.extract_hessian_component(D2psi, H_components, H_components))};
}

template <int dim>
void
run()
{
  using namespace Physics::Invariants;
  namespace TestInvariants = Coupled_Isotropic;

  SymmetricTensor<2, dim> C = unit_symmetric_tensor<dim>();
  for (unsigned int e = 0; e < C.n_independent_components; ++e)
    C.access_raw_entry(e) += (1 + e) * 0.1;
  const SymmetricTensor<2, dim> C_inv = invert(C);

  Tensor<1, dim> H;
  for (unsigned int i = 0; i < H.n_independent_components; ++i)
    H[i] = 10 * (1 + i);

  const auto test_invariant = [&C,
                               &C_inv,
                               &H](const enum InvariantList &    i,
                                   const psi_function_type<dim> &get_psi_ad) {
    const Values<dim> values = compute_derivatives_using_AD(C, H, get_psi_ad);

    std::cout << "\nCmpd Val: " << TestInvariants::Ii(i, C, C_inv, H)
              << "\nExpt Val: " << values.Psi << std::endl;
    std::cout << "\nCmpd dVal: " << TestInvariants::dIi_dC(i, C, C_inv, H)
              << "\nExpt dVal: " << values.dPsi_dC << std::endl;
    std::cout << "\nCmpd dVal: " << TestInvariants::dIi_dH(i, C, C_inv, H)
              << "\nExpt dVal: " << values.dPsi_dH << std::endl;
    std::cout << "\nCmpd ddVal: " << TestInvariants::d2Ii_dC_dC(i, C, C_inv, H)
              << "\nExpt ddVal: " << values.d2Psi_dC_dC << std::endl;
    std::cout << "\nCmpd ddVal: " << TestInvariants::d2Ii_dH_dH(i, C, C_inv, H)
              << "\nExpt ddVal: " << values.d2Psi_dH_dH << std::endl;
    std::cout << "\nCmpd ddVal: " << TestInvariants::d2Ii_dC_dH(i, C, C_inv, H)
              << "\nExpt ddVal: " << values.d2Psi_dC_dH << std::endl;
    std::cout << "\nCmpd ddVal: " << TestInvariants::d2Ii_dH_dC(i, C, C_inv, H)
              << "\nExpt ddVal: " << values.d2Psi_dH_dC << std::endl;

    const double tol = (dim == 2 ? 1e-12 : 1e-9);
    Assert(std::abs(TestInvariants::Ii(i, C, C_inv, H) - values.Psi) < tol,
           ExcMessage("No match in value."));
    Assert((TestInvariants::dIi_dC(i, C, C_inv, H) - values.dPsi_dC).norm() <
             tol,
           ExcMessage("No match in first derivative."));
    Assert((TestInvariants::dIi_dH(i, C, C_inv, H) - values.dPsi_dH).norm() <
             tol,
           ExcMessage("No match in first derivative."));
    Assert((TestInvariants::d2Ii_dC_dC(i, C, C_inv, H) - values.d2Psi_dC_dC)
               .norm() < tol,
           ExcMessage("No match in second derivative."));
    Assert((TestInvariants::d2Ii_dH_dH(i, C, C_inv, H) - values.d2Psi_dH_dH)
               .norm() < tol,
           ExcMessage("No match in second derivative."));
    Assert((TestInvariants::d2Ii_dC_dH(i, C, C_inv, H) - values.d2Psi_dC_dH)
               .norm() < tol,
           ExcMessage("No match in second derivative."));
    Assert((TestInvariants::d2Ii_dH_dC(i, C, C_inv, H) - values.d2Psi_dH_dC)
               .norm() < tol,
           ExcMessage("No match in second derivative."));
  };

  // First invariant
  {
    std::cout << "Checking I1" << std::endl;
    const enum InvariantList invariant = I1;

    const psi_function_type<dim> get_psi_ad =
      [invariant](const SymmetricTensor<2, dim, ADNumberType> &C_ad,
                  const Tensor<1, dim, ADNumberType> &H_ad) -> ADNumberType {
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      const ADNumberType                          psi_ad =
        TestInvariants::Ii(invariant, C_ad, C_inv_ad, H_ad);
      return psi_ad;
    };

    test_invariant(invariant, get_psi_ad);
  }

  // Second invariant
  {
    std::cout << "Checking I2" << std::endl;
    const enum InvariantList invariant = I2;

    const psi_function_type<dim> get_psi_ad =
      [invariant](const SymmetricTensor<2, dim, ADNumberType> &C_ad,
                  const Tensor<1, dim, ADNumberType> &H_ad) -> ADNumberType {
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      const ADNumberType                          psi_ad =
        TestInvariants::Ii(invariant, C_ad, C_inv_ad, H_ad);
      return psi_ad;
    };

    test_invariant(invariant, get_psi_ad);
  }

  // Third invariant
  {
    std::cout << "Checking I3" << std::endl;
    const enum InvariantList invariant = I3;

    const psi_function_type<dim> get_psi_ad =
      [invariant](const SymmetricTensor<2, dim, ADNumberType> &C_ad,
                  const Tensor<1, dim, ADNumberType> &H_ad) -> ADNumberType {
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      const ADNumberType                          psi_ad =
        TestInvariants::Ii(invariant, C_ad, C_inv_ad, H_ad);
      return psi_ad;
    };

    test_invariant(invariant, get_psi_ad);
  }

  // Pseudo third invariant
  {
    std::cout << "Checking pI3" << std::endl;
    const enum InvariantList invariant = pI3;

    const psi_function_type<dim> get_psi_ad =
      [invariant](const SymmetricTensor<2, dim, ADNumberType> &C_ad,
                  const Tensor<1, dim, ADNumberType> &H_ad) -> ADNumberType {
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      const ADNumberType                          psi_ad =
        TestInvariants::Ii(invariant, C_ad, C_inv_ad, H_ad);
      return psi_ad;
    };

    test_invariant(invariant, get_psi_ad);
  }

  // Ninth invariant
  {
    std::cout << "Checking I9" << std::endl;
    const enum InvariantList invariant = I9;

    const psi_function_type<dim> get_psi_ad =
      [invariant](const SymmetricTensor<2, dim, ADNumberType> &C_ad,
                  const Tensor<1, dim, ADNumberType> &H_ad) -> ADNumberType {
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      const ADNumberType                          psi_ad =
        TestInvariants::Ii(invariant, C_ad, C_inv_ad, H_ad);
      return psi_ad;
    };

    test_invariant(invariant, get_psi_ad);
  }

  // Tenth invariant
  {
    std::cout << "Checking I10" << std::endl;
    const enum InvariantList invariant = I10;

    const psi_function_type<dim> get_psi_ad =
      [invariant](const SymmetricTensor<2, dim, ADNumberType> &C_ad,
                  const Tensor<1, dim, ADNumberType> &H_ad) -> ADNumberType {
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      const ADNumberType                          psi_ad =
        TestInvariants::Ii(invariant, C_ad, C_inv_ad, H_ad);
      return psi_ad;
    };

    test_invariant(invariant, get_psi_ad);
  }

  // Eleventh invariant
  {
    std::cout << "Checking 11" << std::endl;
    const enum InvariantList invariant = I11;

    const psi_function_type<dim> get_psi_ad =
      [invariant](const SymmetricTensor<2, dim, ADNumberType> &C_ad,
                  const Tensor<1, dim, ADNumberType> &H_ad) -> ADNumberType {
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      const ADNumberType                          psi_ad =
        TestInvariants::Ii(invariant, C_ad, C_inv_ad, H_ad);
      return psi_ad;
    };

    test_invariant(invariant, get_psi_ad);
  }

  // Pseudo eleventh invariant (a)
  {
    std::cout << "Checking pI11a" << std::endl;
    const enum InvariantList invariant = pI11a;

    const psi_function_type<dim> get_psi_ad =
      [invariant](const SymmetricTensor<2, dim, ADNumberType> &C_ad,
                  const Tensor<1, dim, ADNumberType> &H_ad) -> ADNumberType {
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      const ADNumberType                          psi_ad =
        TestInvariants::Ii(invariant, C_ad, C_inv_ad, H_ad);
      return psi_ad;
    };

    test_invariant(invariant, get_psi_ad);
  }

  // Pseudo eleventh invariant (b)
  {
    std::cout << "Checking pI11b" << std::endl;
    const enum InvariantList invariant = pI11b;

    const psi_function_type<dim> get_psi_ad =
      [invariant](const SymmetricTensor<2, dim, ADNumberType> &C_ad,
                  const Tensor<1, dim, ADNumberType> &H_ad) -> ADNumberType {
      const SymmetricTensor<2, dim, ADNumberType> C_inv_ad = invert(C_ad);
      const ADNumberType                          psi_ad =
        TestInvariants::Ii(invariant, C_ad, C_inv_ad, H_ad);
      return psi_ad;
    };

    test_invariant(invariant, get_psi_ad);
  }

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  run<2>();
  run<3>();

  deallog << "OK" << std::endl;
}
