// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Like q_projection_01 but test project_to_all_subfaces
// Test QProjection for QGaussSimplex.


#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/reference_cell.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int)
{
  AssertThrow(false, ExcNotImplemented());
}

template <>
void
test<2>(const unsigned int n_points)
{
  const int dim = 2;

  QGaussSimplex<dim - 1> quad_ref(n_points);

  // test: QProjector::project_to_all_subfaces
  {
    deallog.push("project_to_all_subfaces");

    const auto quad =
      QProjector<dim>::project_to_all_subfaces(ReferenceCells::Triangle,
                                               quad_ref);

    for (unsigned int q = 0; q < quad.size(); ++q)
      {
        deallog << quad.point(q) << ' ';
        deallog << quad.weight(q) << ' ';
        deallog << std::endl;
      }
    deallog << std::endl;

    deallog << std::endl;
    deallog.pop();
  }


  // test unsigned int DataSetDescriptor::subface
  {
    deallog.push("DataSetDescriptor::subface");
    const unsigned int nq = quad_ref.size();
    for (const unsigned int face_no :
         ReferenceCells::get_simplex<dim>().face_indices())
      for (unsigned int o = 0; o < 2; ++o)
        for (unsigned int subface_no = 0; subface_no < 2; ++subface_no)
          {
            unsigned int offset = QProjector<dim>::DataSetDescriptor::subface(
              ReferenceCells::get_simplex<dim>(), face_no, subface_no, o, nq);
            deallog << "Face:" << face_no << " subface:" << subface_no
                    << " orientation:" << o << " offset:" << offset
                    << std::endl;
          }
    deallog.pop();
  }

  {
    const auto quad =
      QProjector<dim>::project_to_all_subfaces(ReferenceCells::Triangle,
                                               quad_ref);

    const auto print = [&](const unsigned int face_no,
                           const unsigned int subface_no,
                           const unsigned int face_orientation) {
      deallog << "face_no=" << face_no << " subface_no=" << subface_no
              << " face_orientation=" << face_orientation << ':' << std::endl;
      for (unsigned int q = 0,
                        i = QProjector<dim>::DataSetDescriptor::subface(
                          ReferenceCells::Triangle,
                          face_no,
                          subface_no,
                          face_orientation,
                          false,
                          false,
                          quad_ref.size());
           q < quad_ref.size();
           ++q, ++i)
        {
          deallog << quad.point(i) << ' ';
          deallog << quad.weight(i) << ' ';
          deallog << std::endl;
        }
      deallog << std::endl;
    };

    // Loop over faces and loop over subfaces
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 2; ++j)
        {
          print(i, j, 0);
          print(i, j, 1);
        }
  }
}

int
main()
{
  initlog();

  {
    deallog.push("2d-1");
    test<2>(1);
    deallog.pop();
  }
  {
    deallog.push("2d-2");
    test<2>(2);
    deallog.pop();
  }
  {
    deallog.push("2d-3");
    test<2>(3);
    deallog.pop();
  }
}
