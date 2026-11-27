// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check that constructing an FEValues object with a finite element defined
// on a non-hypercube reference cell and a tensor-product quadrature formula
// (QGauss) triggers an assertion, whereas the matching quadrature formula
// works.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"


template <int dim>
void
check(const FiniteElement<dim> &fe, const Quadrature<dim> &quadrature)
{
  deallog << fe.get_name() << " with " << quadrature.size()
          << " quadrature points: ";
  try
    {
      FEValues<dim> fe_values(fe, quadrature, update_values);
      deallog << "OK" << std::endl;
    }
  catch (const ExceptionBase &)
    {
      deallog << "assertion triggered" << std::endl;
    }
}



int
main()
{
  initlog();
  deal_II_exceptions::disable_abort_on_exception();

  for (unsigned int degree = 1; degree <= 2; ++degree)
    {
      check(FE_SimplexP<2>(degree), QGauss<2>(degree + 1));
      check(FE_SimplexP<2>(degree), QGaussSimplex<2>(degree + 1));
      check(FE_SimplexP<3>(degree), QGauss<3>(degree + 1));
      check(FE_SimplexP<3>(degree), QGaussSimplex<3>(degree + 1));
    }

  check(FE_WedgeP<3>(1), QGauss<3>(2));
  check(FE_WedgeP<3>(1), QGaussWedge<3>(2));
  check(FE_PyramidP<3>(1), QGauss<3>(2));
  check(FE_PyramidP<3>(1), QGaussPyramid<3>(2));

  // single-point quadrature formulas and quadrature formulas created from a
  // list of points are used to evaluate at arbitrary points and must still
  // work
  check(FE_SimplexP<2>(1), Quadrature<2>(Point<2>(0.25, 0.25)));
  check(FE_SimplexP<2>(1),
        Quadrature<2>(std::vector<Point<2>>{Point<2>(0.25, 0.25),
                                            Point<2>(0.5, 0.25),
                                            Point<2>(0.25, 0.5)}));
}
