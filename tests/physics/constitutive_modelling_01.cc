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


// Test constitutive models: Isotropic

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/physics/constitutive_modelling.h>
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
  std::function<ADNumberType(const SymmetricTensor<2, dim, ADNumberType> &)>;

template <int dim>
struct Values
{
  double                  Psi;
  SymmetricTensor<2, dim> dPsi_dC;
  SymmetricTensor<4, dim> d2Psi_dC_dC;
};

template <int dim>
Values<dim>
compute_derivatives_using_AD(const SymmetricTensor<2, dim> &C,
                             const psi_function_type<dim> & get_psi_ad)
{
  const FEValuesExtractors::SymmetricTensor<2> C_components(0);
  const unsigned int n_independent_variables = C.n_independent_components;
  ADHelper<dim>      ad_helper(n_independent_variables);

  ad_helper.register_independent_variable(C, C_components);
  const SymmetricTensor<2, dim, ADNumberType> C_ad =
    ad_helper.get_sensitive_variables(C_components);

  const ADNumberType psi_ad = get_psi_ad(C_ad);
  ad_helper.register_dependent_variable(psi_ad);

  Vector<double>     Dpsi(ad_helper.n_independent_variables());
  FullMatrix<double> D2psi(ad_helper.n_independent_variables(),
                           ad_helper.n_independent_variables());

  ad_helper.compute_gradient(Dpsi);
  ad_helper.compute_hessian(D2psi);

  return {ad_helper.compute_value(),
          ad_helper.extract_gradient_component(Dpsi, C_components),
          ad_helper.extract_hessian_component(D2psi,
                                              C_components,
                                              C_components)};
}

template <int dim>
void
run()
{
  using namespace Physics::Invariants;
  using namespace Physics::ConstitutiveModelling;
  using TestInvariants   = Isotropic<dim, double>;
  using TestInvariantsAD = Isotropic<dim, ADNumberType>;

  SymmetricTensor<2, dim> C = unit_symmetric_tensor<dim>();
  for (unsigned int e = 0; e < C.n_independent_components; ++e)
    C.access_raw_entry(e) += (1 + e) * 0.1;
  const SymmetricTensor<2, dim> C_inv = invert(C);

  // For the energy function, we just multiply an easy-to-derive
  // multiplicatively decomposed polynomial function of the invariants, Psi =
  // \Pi_{k}^{n invariants} I_{k}^{3} Raising each invariant to the third power
  // ensures that the invariant value remains in the coefficient when twice
  // differentiated.
  const psi_function_type<dim> get_psi_ad =
    [](const SymmetricTensor<2, dim, ADNumberType> &C_ad) -> ADNumberType {
    ADNumberType psi_ad = 1.0;

    for (const auto i : TestInvariantsAD::valid_invariants())
      psi_ad *= std::pow(TestInvariantsAD::Ii(i, C_ad), 3);

    return psi_ad;
  };

  // Now test the constitutive model. Knowing the form of the
  // energy function allows us to derive the necessary coefficients,
  // with first derivatives
  // dPsi_dIi = 3 I_{i}^{2} * \Pi_{k, k != i}^{n invariants} I_{k}^{3}
  // and second derivatives
  // d2Psi_dIi_dIj = 6 I_{i} * \Pi_{k, k != i}^{n inv.} I_{k}^{3}
  //   when i == j and
  // d2Psi_dIi_dIj = 3 I_{i}^{2} * 3 I_{j}^{2}
  //               * \Pi_{k, k != [i,j]}^{n inv.} I_{k}^{3}
  //   when i != j.
  {
    const Values<dim> values = compute_derivatives_using_AD(C, get_psi_ad);

    using ConstitutiveModel =
      UncoupledConstitutiveModel<TestInvariants>;
    ConstitutiveModel cm;

    for (const auto i : TestInvariants::valid_invariants())
      {
        // Construct the first derivative coefficients
        double dPsi_dIi = 1.0;
        for (const auto kk : TestInvariants::valid_invariants())
          {
            if (kk == i)
              dPsi_dIi *=
                3.0 * std::pow(ConstitutiveModel::InvariantsType::Ii(kk, C), 2);
            else
              dPsi_dIi *=
                std::pow(ConstitutiveModel::InvariantsType::Ii(kk, C), 3);
          }
        cm.set_first_derivative_coefficient(i, dPsi_dIi);

        // Construct the second derivative coefficients
        for (const auto j : TestInvariants::valid_invariants())
          {
            double d2Psi_dIi_dIj = 1.0;

            for (const auto kk : TestInvariants::valid_invariants())
              {
                if (i == j)
                  {
                    if (kk == i)
                      d2Psi_dIi_dIj *=
                        6.0 * ConstitutiveModel::InvariantsType::Ii(kk, C);
                    else
                      d2Psi_dIi_dIj *=
                        std::pow(ConstitutiveModel::InvariantsType::Ii(kk, C),
                                 3);
                  }
                else
                  {
                    if (kk == i || kk == j)
                      d2Psi_dIi_dIj *=
                        3.0 *
                        std::pow(ConstitutiveModel::InvariantsType::Ii(kk, C),
                                 2);
                    else
                      d2Psi_dIi_dIj *=
                        std::pow(ConstitutiveModel::InvariantsType::Ii(kk, C),
                                 3);
                  }
              }

            cm.set_second_derivative_coefficient(i, j, d2Psi_dIi_dIj);
          }
      }

    std::cout << "\nCmpd dPsi_dC: " << cm.get_dPsi_dC(C, C_inv)
              << "\nExpt dPsi_dC: " << values.dPsi_dC << std::endl;
    std::cout << "\nCmpd d2Psi_dC_dC: " << cm.get_d2Psi_dC_dC(C, C_inv)
              << "\nExpt d2Psi_dC_dC: " << values.d2Psi_dC_dC << std::endl;

    const double tol = (dim == 2 ? 1e-12 : 1e-9);
    Assert((cm.get_dPsi_dC(C, C_inv) - values.dPsi_dC).norm() < tol,
           ExcMessage("No match in first derivative."));
    Assert((cm.get_d2Psi_dC_dC(C, C_inv) - values.d2Psi_dC_dC).norm() < tol,
           ExcMessage("No match in second derivative."));

    // Check that we clear with no errors.
    cm.clear();
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
