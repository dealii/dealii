// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2023 by the deal.II authors
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


// Test consistency of ReferenceCell::face_to_cell_vertices(),
// ReferenceCell::vertex(), and ReferenceCell:face_vertex_location().

#include <deal.II/grid/reference_cell.h>

#include <random>

#include "../tests.h"


using namespace dealii;

template <int dim>
void
test(const ReferenceCell &reference_cell)
{
  deallog << "ReferenceCell: " << reference_cell.to_string() << std::endl;
  for (unsigned int f = 0; f < reference_cell.n_faces(); ++f)
    for (unsigned int v = 0;
         v < reference_cell.face_reference_cell(f).n_vertices();
         ++v)
      {
        deallog << "Face=" << f << ", vertex=" << v << ", location="
                << reference_cell.face_vertex_location<dim>(f, v) << std::endl;
        Assert(reference_cell.face_vertex_location<dim>(f, v) ==
                 reference_cell.vertex<dim>(
                   reference_cell.face_to_cell_vertices(
                     f, v, ReferenceCell::default_combined_face_orientation())),
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
