// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
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


#include <deal.II/base/polynomials_integrated_legendre_sz.h>

DEAL_II_NAMESPACE_OPEN

IntegratedLegendreSZ::IntegratedLegendreSZ (const unsigned int k)
  :
  Polynomials::Polynomial<double> (get_coefficients(k))
{}



const std::vector<double>
IntegratedLegendreSZ::get_coefficients (const unsigned int k)
{
  std::vector<double> coefficients(k+1);

  // first two polynomials are hard-coded:
  if (k==0)
    {
      coefficients[0] = -1.;
      return coefficients;
    }
  else if (k==1)
    {
      coefficients[0] = 0.;
      coefficients[1] = 1.;
      return coefficients;
    }

  // General formula is:
  // k*L_{k}(x) = (2*k-3)*x*L_{k-1} - (k-3)*L_{k-2}.
  std::vector<double> coefficients_km2 = get_coefficients(k-2);
  std::vector<double> coefficients_km1 = get_coefficients(k-1);

  const double a = 1.0 / k;
  const double b = 2.0*k - 3.0;
  const double c = k - 3.0;

  // To maintain stability, delay the division (multiplication by a) until the end.
  for (unsigned int i=1; i<=k-2; i++)
    {
      coefficients[i] = b*coefficients_km1[i-1] - c*coefficients_km2[i];
    }

  coefficients[0] = -c*coefficients_km2[0];
  coefficients[k]   = b*coefficients_km1[k-1];
  coefficients[k-1] = b*coefficients_km1[k-2];

  for (unsigned int i=0; i<coefficients.size(); i++)
    {
      coefficients[i] *= a;
    }

  return coefficients;
}



std::vector<Polynomials::Polynomial<double> >
IntegratedLegendreSZ::generate_complete_basis (const unsigned int degree)
{
  std::vector<Polynomials::Polynomial<double> > v;
  v.reserve(degree + 1);
  for (unsigned int i=0; i<=degree; ++i)
    {
      v.push_back (IntegratedLegendreSZ(i));
    }
  return v;
}



DEAL_II_NAMESPACE_CLOSE
