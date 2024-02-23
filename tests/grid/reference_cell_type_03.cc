// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ReferenceCell::barycenter()

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  Triangulation<dim> triangulation;
  GridGenerator::reference_cell(triangulation, reference_cell);

  std::unique_ptr<Quadrature<dim>> q;
  if (reference_cell == ReferenceCells::Line)
    q = std::make_unique<QGauss<dim>>(2);
  else if (reference_cell == ReferenceCells::Triangle)
    q = std::make_unique<QGaussSimplex<dim>>(2);
  else if (reference_cell == ReferenceCells::Quadrilateral)
    q = std::make_unique<QGauss<dim>>(2);
  else if (reference_cell == ReferenceCells::Tetrahedron)
    q = std::make_unique<QGaussSimplex<dim>>(2);
  else if (reference_cell == ReferenceCells::Wedge)
    q = std::make_unique<QGaussWedge<dim>>(2);
  else if (reference_cell == ReferenceCells::Pyramid)
    q = std::make_unique<QGaussPyramid<dim>>(2);
  else if (reference_cell == ReferenceCells::Hexahedron)
    q = std::make_unique<QGauss<dim>>(2);

  FE_Nothing<dim> fe(reference_cell);

  // Set up the objects to compute an integral on the reference cell
  FEValues<dim> fe_values(fe, *q, update_JxW_values | update_quadrature_points);
  fe_values.reinit(triangulation.begin_active());

  double     volume = 0;
  Point<dim> barycenter;
  for (unsigned int i = 0; i < q->size(); ++i)
    {
      volume += fe_values.JxW(i);
      barycenter += fe_values.quadrature_point(i) * fe_values.JxW(i);
    }
  barycenter /= volume;

  deallog << "ReferenceCell: " << reference_cell.to_string() << std::endl;
  deallog << "  computed barycenter = " << barycenter << std::endl;
  deallog << "  self-reported barycenter = " << reference_cell.barycenter<dim>()
          << std::endl;

  Assert((barycenter - reference_cell.barycenter<dim>()).norm() <= 1e-12,
         ExcInternalError());
}

int
main()
{
  initlog();

  {
    deallog.push("1D");
    test<1>(ReferenceCells::Line);
    deallog.pop();
  }

  {
    deallog.push("2D");
    test<2>(ReferenceCells::Quadrilateral);
    test<2>(ReferenceCells::Triangle);
    deallog.pop();
  }

  {
    deallog.push("3D");
    test<3>(ReferenceCells::Tetrahedron);
    test<3>(ReferenceCells::Pyramid);
    test<3>(ReferenceCells::Wedge);
    test<3>(ReferenceCells::Hexahedron);
    deallog.pop();
  }
}
