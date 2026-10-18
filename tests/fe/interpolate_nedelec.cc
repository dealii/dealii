// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nedelec.h>

#include "interpolate_common.h"

// FE_Nedelec<dim>::interpolate(...)

template <int dim>
void
check1(const Function<dim> &f, const unsigned int degree)
{
  FE_Nedelec<dim> fe(degree);
  deallog << fe.get_name() << ' ';
  deallog << fe.get_generalized_support_points().size() << ' ';

  std::vector<double> dofs(fe.dofs_per_cell);

  std::vector<Vector<double>> values(fe.get_generalized_support_points().size(),
                                     Vector<double>(dim));
  f.vector_value_list(fe.get_generalized_support_points(), values);
  fe.convert_generalized_support_point_values_to_dof_values(values, dofs);

  const double error = vector_difference(fe, dofs, f, 0);
  AssertThrow(error < 1e-11,
              ExcMessage("Interpolation must reproduce the polynomial"));
  deallog << " vector " << error << std::endl;
}

// A polynomial in the Nedelec space of order k: component b has degree
// zero in x_b, degree k in x_a, and degree k+1 in the remaining coordinate.
// Its tangential trace vanishes on all edges and differs on the two faces
// normal to a. Thus every face must use its own support-point values.
class FacePolynomial : public Function<3>
{
public:
  FacePolynomial(const unsigned int a,
                 const unsigned int b,
                 const unsigned int k)
    : Function<3>(3)
    , a(a)
    , b(b)
    , k(k)
  {}

  void
  vector_value_list(const std::vector<Point<3>> &points,
                    std::vector<Vector<double>> &values) const override
  {
    const unsigned int c = 3 - a - b;
    for (unsigned int i = 0; i < points.size(); ++i)
      {
        values[i]    = 0;
        values[i][b] = (1. + std::pow(points[i][a], k)) *
                       std::pow(points[i][c], k) * (1. - points[i][c]);
      }
  }

private:
  const unsigned int a, b, k;
};

// Interpolation of an interior shape function must return its unit vector.
// This checks every tensor-index permutation in each component block.
void
check_interior_basis(const unsigned int order)
{
  const FE_Nedelec<3>         fe(order);
  const auto                 &points = fe.get_generalized_support_points();
  std::vector<Vector<double>> values(points.size(), Vector<double>(3));
  std::vector<double>         dofs(fe.n_dofs_per_cell());
  for (unsigned int i = fe.get_first_hex_index(); i < fe.n_dofs_per_cell(); ++i)
    {
      for (unsigned int q = 0; q < points.size(); ++q)
        for (unsigned int c = 0; c < 3; ++c)
          values[q][c] = fe.shape_value_component(i, points[q], c);
      fe.convert_generalized_support_point_values_to_dof_values(values, dofs);
      for (unsigned int j = 0; j < dofs.size(); ++j)
        AssertThrow(std::abs(dofs[j] - (i == j ? 1. : 0.)) < 1e-11,
                    ExcMessage(
                      "Interpolation must preserve interior shape functions"));
    }
  deallog << "Interior basis " << fe.get_name() << " OK" << std::endl;
}

int
main()
{
  initlog();

  {
    Q1WedgeFunction<2, 0, 2> w20;
    check1(w20, 0);
    check1(w20, 1);
    check1(w20, 2);
    Q1WedgeFunction<2, 1, 2> w21;
    check1(w21, 2);
    Q1WedgeFunction<2, 2, 2> w22;
    check1(w22, 2);
    Q1WedgeFunction<2, 3, 2> w23;
    check1(w23, 3);

    PolynomialFunction<2, 1, 2> p21;
    check1(p21, 1);
    PolynomialFunction<2, 2, 2> p22;
    check1(p22, 2);
    PolynomialFunction<2, 3, 2> p23;
    check1(p23, 3);
  }

  {
    // Original issue #328: (xyz)^degree in every component.
    Q1WedgeFunction<3, 0, 3> w0;
    check1(w0, 0);
    Q1WedgeFunction<3, 1, 3> w1;
    check1(w1, 1);
    check1(w1, 2);
    Q1WedgeFunction<3, 2, 3> w2;
    check1(w2, 2);
    Q1WedgeFunction<3, 3, 3> w3;
    check1(w3, 3);
  }

  for (unsigned int k = 1; k <= 3; ++k)
    for (unsigned int a = 0; a < 3; ++a)
      for (unsigned int b = 0; b < 3; ++b)
        if (a != b)
          {
            const FacePolynomial polynomial(a, b, k);
            deallog << "Trace " << a << ' ' << b << ' ';
            check1(polynomial, k);
          }
  check_interior_basis(2);
  check_interior_basis(3);
}
