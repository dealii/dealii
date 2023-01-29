// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

// Test BarycentricPolynomial and BarycentricPolynomials.

#include <deal.II/base/point.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/table.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"

using namespace dealii;

int
main()
{
  initlog();

  BarycentricPolynomial<2> bp2({1, 0, 0}, 1.0);
  deallog << bp2 << std::endl;

  // test some basic algebra with barycentric polynomials
  {
    deallog << "1D:" << std::endl;
    const auto bp1_0 = BarycentricPolynomial<1>::monomial(0);
    const auto bp1_1 = BarycentricPolynomial<1>::monomial(1);

    deallog << "bp1_0 = " << bp1_0 << std::endl;
    deallog << "bp1_1 = " << bp1_1 << std::endl;
    deallog << "bp1_0 * 2 * bp1_1 / 2 = " << bp1_0 * 2 * bp1_1 / 2 << std::endl
            << std::endl;
  }

  {
    deallog << std::endl << "2D:" << std::endl;
    const auto bp2_0 = BarycentricPolynomial<2>::monomial(0) * 2;
    deallog << "bp2_0 = " << bp2_0 << std::endl;

    const auto bp2_1 = 3.0 * BarycentricPolynomial<2>::monomial(1);
    deallog << "bp2_1 = " << bp2_1 << std::endl;

    const auto bp2_2 = BarycentricPolynomial<2>::monomial(2);
    deallog << "bp2_2 = " << bp2_2 << std::endl;

    const auto prod1 = bp2_0 + bp2_1;
    deallog << "bp2_0 + bp2_1 = " << prod1 << std::endl;

    const auto prod2 = prod1 * bp2_0;
    deallog << "(bp2_0 + bp2_1) * bp2_0 = " << prod2 << std::endl;
    deallog << "bp2_0 * bp2_0 + bp2_1 * bp2_0 = "
            << bp2_0 * bp2_0 + bp2_1 * bp2_0 << std::endl;
    deallog << "bp2_1 * bp2_0 + bp2_0 * bp2_0 = "
            << bp2_1 * bp2_0 + bp2_0 * bp2_0 << std::endl;

    // test derivatives
    deallog << "d/dx bp2_0 = " << bp2_0.derivative(0) << std::endl;
    deallog << "d/dy bp2_0 = " << bp2_0.derivative(1) << std::endl;

    deallog << "d/dx bp2_2 = " << bp2_2.derivative(0) << std::endl;
    deallog << "d/dy bp2_2 = " << bp2_2.derivative(1) << std::endl;
  }

  // test various finite element spaces
  {
    deallog << std::endl << "Test with TRI6" << std::endl;

    const auto t1 = BarycentricPolynomial<2>::monomial(0);
    const auto t2 = BarycentricPolynomial<2>::monomial(1);
    const auto t3 = BarycentricPolynomial<2>::monomial(2);

    std::vector<BarycentricPolynomial<2>> p2;
    p2.push_back(t1 * (2 * t1 - 1));
    p2.push_back(t2 * (2 * t2 - 1));
    p2.push_back(t3 * (2 * t3 - 1));
    p2.push_back(4 * t2 * t1);
    p2.push_back(4 * t2 * t3);
    p2.push_back(4 * t3 * t1);

    FE_SimplexP<2> fe(2);
    for (unsigned int i = 0; i < 6; ++i)
      {
        deallog << "p = " << p2[i] << std::endl;
        deallog << "p_x = " << p2[i].derivative(0) << std::endl;
        deallog << "p_y = " << p2[i].derivative(1) << std::endl;
        for (unsigned int j = 0; j < 6; ++j)
          {
            Assert(std::abs(p2[i].value(fe.get_unit_support_points()[j]) -
                            double(i == j)) < 1e-12,
                   ExcInternalError());
          }
        deallog << std::endl;
      }
    deallog << "Test with TRI6 - Success" << std::endl;
  }

  {
    deallog << "Test with TRI10" << std::endl;
    const auto tri10 = BarycentricPolynomials<2>::get_fe_p_basis(3);

    FE_SimplexP<2> fe(3);
    const auto &   points = fe.get_unit_support_points();
    for (unsigned int i = 0; i < 10; ++i)
      {
        Assert(points.size() == 10, ExcInternalError());
        for (unsigned int j = 0; j < 10; ++j)
          {
            Assert(std::abs(tri10.compute_value(i, points[j]) -
                            double(i == j)) < 1e-12,
                   ExcInternalError());

            // third derivatives should be constant
            Assert((tri10.compute_3rd_derivative(i, points[0]) -
                    tri10.compute_3rd_derivative(i, points[j]))
                       .norm() == 0.0,
                   ExcInternalError());

            Assert(tri10.compute_4th_derivative(i, points[j]).norm() == 0.0,
                   ExcInternalError());
          }
      }
    deallog << "Test with TRI10 - Success" << std::endl;
  }

  {
    deallog << "Test with TET4" << std::endl;
    const auto tet4 = BarycentricPolynomials<3>::get_fe_p_basis(1);

    FE_SimplexP<3> fe(1);
    const auto &   points = fe.get_unit_support_points();
    for (unsigned int i = 0; i < 4; ++i)
      {
        Assert(points.size() == 4, ExcInternalError());
        for (unsigned int j = 0; j < 4; ++j)
          {
            Assert(std::abs(tet4.compute_value(i, points[j]) - double(i == j)) <
                     1e-12,
                   ExcInternalError());

            // first derivatives should be constant
            Assert((tet4.compute_grad(i, points[0]) -
                    tet4.compute_grad(i, points[j]))
                       .norm() == 0.0,
                   ExcInternalError());
            Assert(tet4.compute_2nd_derivative(i, points[j]).norm() == 0.0,
                   ExcInternalError());
            Assert(tet4.compute_3rd_derivative(i, points[j]).norm() == 0.0,
                   ExcInternalError());
            Assert(tet4.compute_4th_derivative(i, points[j]).norm() == 0.0,
                   ExcInternalError());
          }
      }
    deallog << "Test with TET4 - Success" << std::endl;
  }

  {
    deallog << "Test with TET10" << std::endl;
    const auto tet10 = BarycentricPolynomials<3>::get_fe_p_basis(2);

    FE_SimplexP<3> fe(2);
    const auto &   points = fe.get_unit_support_points();
    for (unsigned int i = 0; i < 10; ++i)
      {
        Assert(points.size() == 10, ExcInternalError());
        for (unsigned int j = 0; j < 10; ++j)
          {
            Assert(std::abs(tet10.compute_value(i, points[j]) -
                            double(i == j)) < 1e-12,
                   ExcInternalError());

            // second derivatives should be constant
            Assert((tet10.compute_2nd_derivative(i, points[0]) -
                    tet10.compute_2nd_derivative(i, points[j]))
                       .norm() == 0.0,
                   ExcInternalError());

            Assert(tet10.compute_3rd_derivative(i, points[j]).norm() == 0.0,
                   ExcInternalError());
            Assert(tet10.compute_4th_derivative(i, points[j]).norm() == 0.0,
                   ExcInternalError());
          }
      }
    deallog << "Test with TET10 - Success" << std::endl;
  }

  {
    deallog << "Test with TET20" << std::endl;
    const auto tet20 = BarycentricPolynomials<3>::get_fe_p_basis(3);

    FE_SimplexP<3> fe(3);
    const auto &   points = fe.get_unit_support_points();
    for (unsigned int i = 0; i < 20; ++i)
      {
        Assert(points.size() == 20, ExcInternalError());
        for (unsigned int j = 0; j < 20; ++j)
          {
            Assert(std::abs(tet20.compute_value(i, points[j]) -
                            double(i == j)) < 1e-12,
                   ExcInternalError());

            // third derivatives should be constant
            Assert((tet20.compute_3rd_derivative(i, points[0]) -
                    tet20.compute_3rd_derivative(i, points[j]))
                       .norm() == 0.0,
                   ExcInternalError());

            Assert(tet20.compute_4th_derivative(i, points[j]).norm() == 0.0,
                   ExcInternalError());
          }
      }
    deallog << "Test with TET20 - Success" << std::endl;
  }
}
