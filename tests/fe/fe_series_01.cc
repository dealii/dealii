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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// test Fourier expansion in 1D for some simple functions.
// Below is the MWE in Maxima
/*********************************************************
a:0;
b:1;
nmax:3;
define(Phi(n), exp(-2*%i*%pi*n*x/(b-a)));
f:x;
C0:integrate(f*conjugate(Phi(0)),x,a,b)/(b-a);
define(C(n),integrate(f*conjugate(Phi(n)),x,a,b)/(b-a));
fullratsimp(map(C,makelist(i,i,1,nmax)));
fs(nmax):=C0+sum(realpart(conjugate(C(m))*Phi(-m)+C(m)*Phi(m)),m,1,nmax);
plot2d([f,fs(0),fs(1),fs(2),fs(3)],[x,0,1]);
*********************************************************/

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <iostream>

#include "../tests.h"


void
test_1d()
{
  const unsigned int N = 4;
  // exact values obtained by Maxima:
  std::vector<std::complex<double>> exact(N);
  const std::complex<double>        I(0, 1.);
  exact[0] = std::complex<double>(1. / 2, 0.);
  exact[1] = -I / (2 * numbers::PI);
  exact[2] = -I / (4 * numbers::PI);
  exact[3] = -I / (6 * numbers::PI);
  //
  const unsigned int    dim = 1;
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim>  q_collection;

  // linear FE
  fe_collection.push_back(FE_Q<dim>(1));

  QGauss<1>      base_quadrature(6);
  QIterated<dim> quadrature(base_quadrature, N);
  q_collection.push_back(quadrature);

  const std::vector<unsigned int> n_coefficients_per_direction(1, N);

  FESeries::Fourier<dim> fourier(n_coefficients_per_direction,
                                 fe_collection,
                                 q_collection);

  Vector<double> local_dof_values(2);
  local_dof_values[0]                     = 0;
  local_dof_values[1]                     = 1.;
  const unsigned int cell_active_fe_index = 0;

  Table<dim, std::complex<double>> fourier_coefficients;
  fourier_coefficients.reinit(TableIndices<1>(N));

  fourier.calculate(local_dof_values,
                    cell_active_fe_index,
                    fourier_coefficients);

  deallog << "calculated:" << std::endl;
  for (unsigned int i = 0; i < N; i++)
    deallog << fourier_coefficients[i].real() << " "
            << fourier_coefficients[i].imag() << std::endl;
  deallog << "exact:" << std::endl;
  for (unsigned int i = 0; i < N; i++)
    deallog << exact[i].real() << " " << exact[i].imag() << std::endl;
}



int
main()
{
  initlog();

  test_1d();
}
