// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// Check the mass matrix of the HermiteLikeInterpolation

#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"


int
main()
{
  initlog();

  for (unsigned int degree = 1; degree < 11; ++degree)
    {
      auto poly =
        Polynomials::HermiteLikeInterpolation::generate_complete_basis(degree);
      QGauss<1>          quad(degree + 1);
      FullMatrix<double> mat(degree + 1, degree + 1);
      for (unsigned int i = 0; i < degree + 1; ++i)
        for (unsigned int j = 0; j < degree + 1; ++j)
          {
            double sum = 0;
            for (unsigned int q = 0; q < quad.size(); ++q)
              sum += poly[i].value(quad.point(q)[0]) *
                     poly[j].value(quad.point(q)[0]) * quad.weight(q);
            if (std::abs(sum) > 1e-15)
              mat(i, j) = sum;
          }

      deallog << "Mass matrix degree " << degree << std::endl;
      mat.print_formatted(deallog.get_file_stream());
      deallog << std::endl;
    }
}
