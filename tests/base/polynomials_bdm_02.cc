// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Tests 3D polynomials_BDM grad grad at a collection of points on the unit
// square

#include <deal.II/base/job_identifier.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <vector>

#include "../tests.h"

template <int dim>
void
plot(const PolynomialsBDM<dim> &poly)
{
  const PolynomialSpace<dim> legendre_poly_space =
    Polynomials::Legendre::generate_complete_basis(poly.degree() - 1);

  const Point<3> p0(0, 0, 0);
  const Point<3> p1(0.25, 0.5, 0.75);
  const Point<3> p2(0.75, 0.5, 0.25);
  const Point<3> p3(0.75, 0.25, 0.5);
  const Point<3> p4(1.0, 1.0, 1.0);

  std::vector<Point<3>> points;
  points.push_back(p0);
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);
  points.push_back(p4);

  std::vector<Tensor<1, dim>> values;
  std::vector<Tensor<2, dim>> grads;
  std::vector<Tensor<3, dim>> grads2(poly.n());
  std::vector<Tensor<4, dim>> thirds;
  std::vector<Tensor<5, dim>> fourths;

  const unsigned int n_sub = legendre_poly_space.n();

  for (unsigned int k = 0; k < points.size(); ++k)
    {
      if (k % (poly.degree() + 3) == 0)
        deallog << "BDM" << poly.degree() - 1 << '<' << dim << '>' << std::endl;

      unsigned int start = dim * n_sub;

      deallog << "BDM" << poly.degree() - 1 << '<' << dim << '>' << points[k]
              << std::endl;
      poly.evaluate(points[k], values, grads, grads2, thirds, fourths);

      for (unsigned int i = 0; i < poly.degree(); ++i, start += dim)
        for (unsigned int j = 0; j < dim; ++j)
          {
            for (unsigned int d1 = 0; d1 < dim; ++d1)
              for (unsigned int d2 = 0; d2 < dim; ++d2)
                for (unsigned int d3 = 0; d3 < dim; ++d3)
                  deallog << '\t' << grads2[start + j][d1][d2][d3];
          }
      deallog << std::endl;
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  PolynomialsBDM<3> p31(1);
  PolynomialsBDM<3> p32(2);
  PolynomialsBDM<3> p33(3);

  plot(p31);
  plot(p32);
  plot(p33);
}
