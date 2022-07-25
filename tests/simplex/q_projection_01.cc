// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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


// Test QProjection for QGaussSimplex.


#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test(const unsigned int)
{
  AssertThrow(false, ExcNotImplemented())
}

template <>
void
test<2>(const unsigned int n_points)
{
  const int dim = 2;

  QGaussSimplex<dim - 1> quad_ref(n_points);

  const auto quad =
    QProjector<dim>::project_to_all_faces(ReferenceCells::Triangle, quad_ref);

  const auto print = [&](const unsigned int face_no,
                         const bool         face_orientation) {
    deallog << "face_no=" << face_no
            << " face_orientation=" << (face_orientation ? "true" : "false")
            << ':' << std::endl;
    for (unsigned int
           q = 0,
           i =
             QProjector<dim>::DataSetDescriptor::face(ReferenceCells::Triangle,
                                                      face_no,
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

  for (unsigned int i = 0; i < 3; ++i)
    {
      print(i, true);  // 1
      print(i, false); // 0
    }
}

template <>
void
test<3>(const unsigned int n_points)
{
  const int dim = 3;

  QGaussSimplex<dim - 1> quad_ref(n_points);

  const auto quad =
    QProjector<dim>::project_to_all_faces(ReferenceCells::Tetrahedron,
                                          quad_ref);

  const auto print = [&](const unsigned int face_no,
                         const bool         face_orientation,
                         const bool         face_flip,
                         const bool         face_rotation) {
    deallog << "face_no=" << face_no
            << " face_orientation=" << (face_orientation ? "true" : "false")
            << " face_flip=" << (face_flip ? "true" : "false")
            << " face_rotation=" << (face_rotation ? "true" : "false") << ':'
            << std::endl;
    for (unsigned int q = 0,
                      i = QProjector<dim>::DataSetDescriptor::face(
                        ReferenceCells::Tetrahedron,
                        face_no,
                        face_orientation,
                        face_flip,
                        face_rotation,
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

  for (unsigned int i = 0; i < 4; ++i)
    {
      print(i, true, false, false);  // 1
      print(i, true, true, false);   // 3
      print(i, true, false, true);   // 5
      print(i, false, false, false); // 0
      print(i, false, true, false);  // 2
      print(i, false, false, true);  // 4
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

  {
    deallog.push("3d-1");
    test<3>(1);
    deallog.pop();
  }
  {
    deallog.push("3d-4");
    test<3>(2);
    deallog.pop();
  }
  {
    deallog.push("3d-10");
    test<3>(3);
    deallog.pop();
  }
}
