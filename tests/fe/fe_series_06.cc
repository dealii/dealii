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


// test Fourier expansion in 2D for a given vector of local DoF values.
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <iostream>

#include "../tests.h"


void
test_2d()
{
  const unsigned int    dim = 2;
  const unsigned int    N   = 7;
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim>  fourier_q_collection;

  for (unsigned int degree = 2; degree <= N; ++degree)
    {
      fe_collection.push_back(FE_Q<dim>(degree));
    }

  QGauss<1>      base_quadrature(2);
  QIterated<dim> quadrature(base_quadrature, N);
  for (unsigned int i = 0; i < fe_collection.size(); i++)
    fourier_q_collection.push_back(quadrature);

  FESeries::Fourier<dim> fourier(N, fe_collection, fourier_q_collection);
  Table<dim, std::complex<double>> fourier_coefficients;
  fourier_coefficients.reinit(N, N);

  Vector<double> local_dof_values(9);
  double         dofs[] = {0.0000000000000000e+00,
                   0.0000000000000000e+00,
                   0.0000000000000000e+00,
                   2.3801522930483391e-04,
                   0.0000000000000000e+00,
                   1.1949018981806140e-04,
                   0.0000000000000000e+00,
                   1.1949019042912971e-04,
                   5.9982796422221083e-05};
  for (unsigned int i = 0; i < 9; i++)
    local_dof_values[i] = dofs[i];

  const unsigned int cell_active_fe_index = 0;
  fourier.calculate(local_dof_values,
                    cell_active_fe_index,
                    fourier_coefficients);

  for (unsigned int i = 0; i < fourier_coefficients.size(0); i++)
    for (unsigned int j = 0; j < fourier_coefficients.size(1); j++)
      if ((i * i + j * j < N * N) && (i * i + j * j > 0))
        deallog << (i * i + j * j) << " : " << fourier_coefficients(i, j)
                << std::endl;
}



int
main()
{
  initlog();

  test_2d();
}
