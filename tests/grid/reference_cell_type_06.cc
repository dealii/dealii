// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test ReferenceCell::face_measure()

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/q_collection.h>

#include "../tests.h"


template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  Triangulation<dim> triangulation;
  GridGenerator::reference_cell(triangulation, reference_cell);

  FE_Nothing<dim> fe(reference_cell);

  hp::QCollection<dim - 1> quadratures;
  // TODO: MappingQ asserts that face quadrature collections have exactly one
  // entry. For now just special-case those.
  if (reference_cell.is_hyper_cube())
    quadratures.push_back(reference_cell.face_reference_cell(0)
                            .template get_gauss_type_quadrature<dim - 1>(2));
  else
    for (unsigned int face_no = 0; face_no < reference_cell.n_faces();
         ++face_no)
      quadratures.push_back(reference_cell.face_reference_cell(face_no)
                              .template get_gauss_type_quadrature<dim - 1>(2));

  // Set up the objects to compute an integral on the reference cell
  FEFaceValues<dim> fe_face_values(fe, quadratures, update_JxW_values);

  for (unsigned int face_no = 0; face_no < reference_cell.n_faces(); ++face_no)
    {
      fe_face_values.reinit(triangulation.begin_active(), face_no);

      double measure = 0;
      for (unsigned int i = 0; i < fe_face_values.n_quadrature_points; ++i)
        measure += fe_face_values.JxW(i);

      deallog << "ReferenceCell: " << reference_cell.to_string() << std::endl;
      deallog << "  computed face measure = " << measure << std::endl;
      deallog << "  self-reported face measure = "
              << reference_cell.face_measure(face_no) << std::endl;

      Assert(std::fabs(measure -
                       triangulation.begin_active()->face(face_no)->measure()) <
               1e-12,
             ExcInternalError());
      Assert(std::fabs(measure - reference_cell.face_measure(face_no)) < 1e-12,
             ExcInternalError());
    }
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
