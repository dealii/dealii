// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test PolynomialsWedge by comparing it to the old implementation


#include <deal.II/base/point.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/polynomials_wedge.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"

class PolynomialWedge
{
public:
  PolynomialWedge(const unsigned int degree_in)
    : degree(degree_in)
    , poly_tri(BarycentricPolynomials<2>::get_fe_p_basis(degree_in))
    , poly_line(BarycentricPolynomials<1>::get_fe_p_basis(degree_in)){};


  double
  compute_value(const unsigned int i, const Point<3> &p) const
  {
    const auto pair = degree == 1 ? wedge_table_1[i] : wedge_table_2[i];

    const Point<2> p_tri(p[0], p[1]);
    const auto     v_tri = poly_tri.compute_value(pair[0], p_tri);

    const Point<1> p_line(p[2]);
    const auto     v_line = poly_line.compute_value(pair[1], p_line);

    return v_tri * v_line;
  };

  Tensor<1, 3>
  compute_grad(const unsigned int i, const Point<3> &p) const
  {
    const auto pair = degree == 1 ? wedge_table_1[i] : wedge_table_2[i];

    const Point<2> p_tri(p[0], p[1]);
    const auto     v_tri = poly_tri.compute_value(pair[0], p_tri);
    const auto     g_tri = poly_tri.compute_grad(pair[0], p_tri);

    const Point<1> p_line(p[2]);
    const auto     v_line = poly_line.compute_value(pair[1], p_line);
    const auto     g_line = poly_line.compute_grad(pair[1], p_line);

    Tensor<1, 3> grad;
    grad[0] = g_tri[0] * v_line;
    grad[1] = g_tri[1] * v_line;
    grad[2] = v_tri * g_line[0];

    return grad;
  };

private:
  const unsigned int degree;

  const BarycentricPolynomials<2> poly_tri;
  const BarycentricPolynomials<1> poly_line;


  const dealii::ndarray<unsigned int, 6, 2> wedge_table_1{
    {{{0, 0}}, {{1, 0}}, {{2, 0}}, {{0, 1}}, {{1, 1}}, {{2, 1}}}};

  const dealii::ndarray<unsigned int, 18, 2> wedge_table_2{{{{0, 0}},
                                                            {{1, 0}},
                                                            {{2, 0}},
                                                            {{0, 1}},
                                                            {{1, 1}},
                                                            {{2, 1}},
                                                            {{3, 0}},
                                                            {{4, 0}},
                                                            {{5, 0}},
                                                            {{3, 1}},
                                                            {{4, 1}},
                                                            {{5, 1}},
                                                            {{0, 2}},
                                                            {{1, 2}},
                                                            {{2, 2}},
                                                            {{3, 2}},
                                                            {{4, 2}},
                                                            {{5, 2}}}};
};

template <int dim>
void
test(const unsigned int degree)
{
  FE_WedgeP<dim> fe_wedge(degree);
  const auto     support_points = fe_wedge.get_unit_support_points();
  deallog << "N support points: " << support_points.size() << " at degree "
          << degree << std::endl;
  const auto poly = ScalarLagrangePolynomialWedge(degree, support_points);

  const auto poly_legacy_constructor =
    ScalarLagrangePolynomialWedge<dim>(degree);

  const PolynomialWedge poly_wedge(degree);

  std::vector<Point<dim>> points;
  for (const auto &p : support_points)
    points.emplace_back(p);
  // Add some random points
  points.emplace_back(Point<dim>(0.36658, 0.58775, 0.21455));
  points.emplace_back(Point<dim>(0.64464, 0.3546, 0.3246));
  points.emplace_back(Point<dim>(0.0, 0.0, 0.999));

  // Add points from the quadrature
  QGaussWedge<dim> quad(2);
  for (unsigned int i = 0; i < quad.size(); ++i)
    points.emplace_back(quad.point(i));

  for (unsigned int i = 0; i < points.size(); ++i)
    for (unsigned int j = 0; j < fe_wedge.n_dofs_per_cell(); ++j)
      {
        const auto v1 = poly.compute_value(j, points[i]);
        const auto v2 = poly_wedge.compute_value(j, points[i]);
        const auto v3 = poly_legacy_constructor.compute_value(j, points[i]);
        if (std::abs(v2 - v1) < 1e-9)
          deallog << "ok ";
        else
          deallog << "Failure by value!!! " << i << " " << j << " " << v2 - v1
                  << std::endl;
        if (std::abs(v3 - v1) < 1e-9)
          deallog << "ok ";
        else
          deallog << "Failure by value (legacy constructor)!!! " << i << " "
                  << j << " " << v3 << " " << v1 << std::endl;


        const auto g1 = poly.compute_grad(j, points[i]);
        const auto g2 = poly_wedge.compute_grad(j, points[i]);
        const auto g3 = poly_legacy_constructor.compute_grad(j, points[i]);
        for (unsigned int d = 0; d < dim; ++d)
          {
            if (std::abs((g2 - g1)[d]) < 1e-9)
              deallog << "ok ";
            else
              deallog << "Failure in gradient!!! " << i << " " << j << " "
                      << g2[d] << " " << g1[d] << " " << d << std::endl;
          }
        for (unsigned int d = 0; d < dim; ++d)
          {
            if (std::abs((g3 - g1)[d]) < 1e-9)
              deallog << "ok ";
            else
              deallog << "Failure in gradient (legacy constructor)!!! " << i
                      << " " << j << " " << g3[d] << " " << g1[d] << " " << d
                      << std::endl;
          }
        deallog << std::endl;
      }
  deallog << std::endl;
}


int
main()
{
  initlog();
  {
    deallog.push("3d-1");
    test<3>(1);
    deallog.pop();

    deallog.push("3d-2");
    test<3>(2);
    deallog.pop();
  }
}
