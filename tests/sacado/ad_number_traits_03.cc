// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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


// Test that all of the nested Sacado number extractors return the correct
// information

#include <deal.II/differentiation/ad/ad_number_types.h>
#include <deal.II/differentiation/ad/sacado_number_types.h>

#include <complex>
#include <fstream>
#include <iomanip>
#include <type_traits>

#include "../tests.h"

namespace AD = Differentiation::AD;

template <typename Number>
Number
func(const Number &a, const Number &b, const Number &c)
{
  const Number r = c * std::log(b + 1.) / std::sin(a);
  return r;
}

template <typename Number>
void
print(const Number &r,
      const Number &drda,
      const Number &drdb,
      const Number &d2rda2,
      const Number &d2rdb2,
      const Number &d2rdadb,
      const Number &d2rdbda)
{
  deallog << "        r = " << r << std::endl
          << "    dr/da = " << drda << std::endl
          << "    dr/db = " << drdb << std::endl
          << "d^2r/da^2 = " << d2rda2 << std::endl
          << "d^2r/db^2 = " << d2rdb2 << std::endl
          << "d^2r/dadb = " << d2rdadb << std::endl
          << "d^2r/dbda = " << d2rdbda << std::endl;
}

int
main()
{
  initlog();

  // Values of function arguments
  const double a = M_PI / 4;
  const double b = 2.0;
  const double c = 3.0;

  // Number of independent variables
  const int num_deriv = 2;

  // Fad objects
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> afad(num_deriv, 0, a);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> bfad(num_deriv, 1, b);
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> cfad = c;

  afad.val() = Sacado::Fad::DFad<double>(num_deriv, 0, a);
  bfad.val() = Sacado::Fad::DFad<double>(num_deriv, 1, b);

  // AD typedefs
  typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double>> ADNumberType;
  typedef typename ADNumberType::value_type            ADDerivativeType;
  typedef typename ADNumberType::scalar_type ADScalarType; // == double

  // Compute function and derivative with AD
  const ADNumberType rfad = func(afad, bfad, cfad);

  deallog.push("Native");
  {
    // Extract value and derivatives
    const double r       = rfad.val().val(); // r
    const double drda    = rfad.dx(0).val(); // dr/da
    const double drdb    = rfad.dx(1).val(); // dr/db
    const double d2rda2  = rfad.dx(0).dx(0); // d^2r/da^2
    const double d2rdadb = rfad.dx(0).dx(1); // d^2r/dadb
    const double d2rdbda = rfad.dx(1).dx(0); // d^2r/dbda
    const double d2rdb2  = rfad.dx(1).dx(1); // d^2/db^2

    print(r, drda, drdb, d2rda2, d2rdb2, d2rdadb, d2rdbda);
  }
  deallog.pop();

  deallog.push("Extracted");
  {
    // Extract value and derivatives
    const ADScalarType r =
      AD::internal::ExtractData<ADNumberType>::value(rfad); // r
    const ADDerivativeType drda_ad =
      AD::internal::ExtractData<ADNumberType>::directional_derivative(rfad, 0);
    const ADDerivativeType drdb_ad =
      AD::internal::ExtractData<ADNumberType>::directional_derivative(rfad, 1);
    const ADScalarType drda =
      AD::internal::ExtractData<ADDerivativeType>::value(drda_ad); // dr/da
    const ADScalarType drdb =
      AD::internal::ExtractData<ADDerivativeType>::value(drdb_ad); // dr/db
    const ADScalarType d2rda2 =
      AD::internal::ExtractData<ADDerivativeType>::directional_derivative(
        drda_ad, 0); // d^2r/da^2
    const ADScalarType d2rdadb =
      AD::internal::ExtractData<ADDerivativeType>::directional_derivative(
        drda_ad, 1); // d^2r/dadb
    const ADScalarType d2rdbda =
      AD::internal::ExtractData<ADDerivativeType>::directional_derivative(
        drdb_ad, 0); // d^2r/dbda
    const ADScalarType d2rdb2 =
      AD::internal::ExtractData<ADDerivativeType>::directional_derivative(
        drdb_ad, 1); // d^2/db^2

    print(r, drda, drdb, d2rda2, d2rdb2, d2rdadb, d2rdbda);
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
