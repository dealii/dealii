// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Verify that the quadrature rules created by QProjector::project_to_face() and
// QProjector::project_to_all_faces() are correct. I calculated the reference
// values by-hand.

#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"

template <int dim, std::size_t N, typename F>
void
test(const ReferenceCell             reference_cell,
     const hp::QCollection<dim - 1> &quadrature,
     const std::array<double, N>    &integrals,
     const F                        &integrand)
{
  Assert(reference_cell.n_faces() == N, ExcInternalError());
  const auto all_quadratures =
    QProjector<dim>::project_to_all_faces(reference_cell, quadrature);
  for (const unsigned int face_no : reference_cell.face_indices())
    for (types::geometric_orientation combined_orientation = 0;
         combined_orientation < reference_cell.n_face_orientations(face_no);
         ++combined_orientation)
      {
        const auto offset = QProjector<dim>::DataSetDescriptor::face(
          reference_cell, face_no, combined_orientation, quadrature);
        const Quadrature<dim> face_quadrature =
          QProjector<dim>::project_to_face(
            reference_cell,
            quadrature[quadrature.size() == 1 ? 0 : face_no],
            face_no,
            combined_orientation);
        double integral = 0;
        for (unsigned int qp_n = 0; qp_n < face_quadrature.size(); ++qp_n)
          {
            integral += integrand(face_quadrature.point(qp_n)) *
                        face_quadrature.weight(qp_n);
            Assert(std::abs((all_quadratures.point(offset + qp_n) -
                             face_quadrature.point(qp_n))
                              .norm()) < 1e-14,
                   ExcInternalError());
            Assert(std::abs(all_quadratures.weight(offset + qp_n) -
                            face_quadrature.weight(qp_n)) < 1e-14,
                   ExcInternalError());
          }
        Assert(std::abs(integral - integrals[face_no]) <
                 std::max(std::abs(integral), 1.0) * 1e-14,
               ExcInternalError());
        deallog << "(face_no, combined_orientation) = (" << face_no << ", "
                << int(combined_orientation) << ")" << std::endl
                << "value = " << integral << std::endl;
      }
}

int
main()
{
  initlog();

  // test that we can correctly integrate (x + 1)^2 on each face
  {
    const QGauss<0> quadrature(1);
    const auto      integrand = [&](const Point<1> &p) {
      const double u = p[0] + 1;
      return u * u;
    };
    deallog.push("Line");
    {
      const std::array<double, 2> integrals{{1.0, 4.0}};
      test<1>(ReferenceCells::Line,
              hp::QCollection<0>(quadrature),
              integrals,
              integrand);
    }
    deallog.pop();
  }

  // test that we can correctly integrate (x + 1)^2 * (y + 2)^2 on each face
  {
    const QGauss<1> quadrature(4);
    const auto      integrand = [&](const Point<2> &p) {
      const double u = p[0] + 1, v = p[1] + 2;
      return u * u * v * v;
    };
    deallog.push("Triangle");
    {
      const std::array<double, 3> integrals{
        {28.0 / 3.0, 203.0 / 15.0 * std::sqrt(2.0), 19.0 / 3.0}};
      test<2>(ReferenceCells::Triangle,
              hp::QCollection<1>(quadrature),
              integrals,
              integrand);
    }
    deallog.pop();

    deallog.push("Quadrilateral");
    {
      const std::array<double, 4> integrals{
        {19.0 / 3.0, 76.0 / 3.0, 28.0 / 3.0, 21.0}};
      test<2>(ReferenceCells::Quadrilateral,
              hp::QCollection<1>(quadrature),
              integrals,
              integrand);
    }
    deallog.pop();
  }

  // test that we can correctly integrate (x + 1)^2 (y + 2)^2 + (z + 3)^2 on
  // each face
  {
    const QGauss<2>        quadrature_quadrilateral(4);
    const QGaussSimplex<2> quadrature_triangle(4);
    const auto             integrand = [&](const Point<3> &p) {
      const double u = p[0] + 1.0, v = p[1] + 2.0, w = p[2] + 3.0;
      return u * u * v * v + w * w;
    };
    deallog.push("Tetrahedron");
    {
      const std::array<double, 4> integrals{{421.0 / 45.0,
                                             37.0 / 4.0,
                                             25.0 / 3.0,
                                             1879.0 / 180.0 * std::sqrt(3.0)}};
      test<3>(ReferenceCells::Tetrahedron,
              hp::QCollection<2>(quadrature_triangle),
              integrals,
              integrand);
    }
    deallog.pop();

    deallog.push("Pyramid");
    {
      const std::array<double, 5> integrals{{532.0 / 9.0,
                                             533.0 / 45.0 * std::sqrt(2.0),
                                             1037.0 / 45.0 * std::sqrt(2.0),
                                             596.0 / 45.0 * std::sqrt(2.0),
                                             884.0 / 45.0 * std::sqrt(2.0)}};
      test<3>(ReferenceCells::Pyramid,
              hp::QCollection<2>(quadrature_quadrilateral,
                                 quadrature_triangle,
                                 quadrature_triangle,
                                 quadrature_triangle,
                                 quadrature_triangle),
              integrals,
              integrand);
    }
    deallog.pop();

    const double X = 0.0;
    deallog.push("Wedge");
    {
      const std::array<double, 5> integrals{{421.0 / 45.0,
                                             1157.0 / 90.0,
                                             65.0 / 3.0,
                                             388.0 / 15.0 * std::sqrt(2.0),
                                             56.0 / 3.0}};
      test<3>(ReferenceCells::Wedge,
              hp::QCollection<2>(quadrature_triangle,
                                 quadrature_triangle,
                                 quadrature_quadrilateral,
                                 quadrature_quadrilateral,
                                 quadrature_quadrilateral),
              integrals,
              integrand);
    }
    deallog.pop();

    deallog.push("Hexahedron");
    {
      const std::array<double, 6> integrals{
        {56 / 3.0, 113.0 / 3.0, 65 / 3.0, 100 / 3.0, 214.0 / 9.0, 277 / 9.0}};
      test<3>(ReferenceCells::Hexahedron,
              hp::QCollection<2>(quadrature_quadrilateral),
              integrals,
              integrand);
    }
    deallog.pop();
  }
}
