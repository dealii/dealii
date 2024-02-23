// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// output all integer values and functions of GeometryInfo

#include <deal.II/base/geometry_info.h>

#include "../tests.h"


template <int dim>
void
test()
{
  deallog << "max_children_per_cell "
          << GeometryInfo<dim>::max_children_per_cell << std::endl;
  deallog << "faces_per_cell    " << GeometryInfo<dim>::faces_per_cell
          << std::endl;
  deallog << "max_children_per_face "
          << GeometryInfo<dim>::max_children_per_face << std::endl;
  deallog << "vertices_per_cell " << GeometryInfo<dim>::vertices_per_cell
          << std::endl;
  deallog << "lines_per_cell    " << GeometryInfo<dim>::lines_per_cell
          << std::endl;
  deallog << "quads_per_cell    " << GeometryInfo<dim>::quads_per_cell
          << std::endl;
  deallog << "hexes_per_cell    " << GeometryInfo<dim>::hexes_per_cell
          << std::endl;

  deallog << "vertices_per_face " << GeometryInfo<dim>::vertices_per_face
          << std::endl;
  deallog << "lines_per_face    " << GeometryInfo<dim>::lines_per_face
          << std::endl;
  deallog << "quads_per_face    " << GeometryInfo<dim>::quads_per_face
          << std::endl;

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    deallog << "face normal" << f << ' '
            << (GeometryInfo<dim>::unit_normal_orientation[f] > 0. ? '+' : '-')
            << 'x' << GeometryInfo<dim>::unit_normal_direction[f] << std::endl;

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    {
      deallog << "face_children" << f << "[true ]";
      for (unsigned int v = 0; v < GeometryInfo<dim>::max_children_per_face;
           ++v)
        deallog << ' '
                << GeometryInfo<dim>::child_cell_on_face(
                     RefinementCase<dim>::isotropic_refinement, f, v, true);
      deallog << std::endl;
      deallog << "face_children" << f << "[false]";
      for (unsigned int v = 0; v < GeometryInfo<dim>::max_children_per_face;
           ++v)
        deallog << ' '
                << GeometryInfo<dim>::child_cell_on_face(
                     RefinementCase<dim>::isotropic_refinement, f, v, false);
      deallog << std::endl;
    }

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    {
      deallog << "face_vertices" << f << "[true ]";
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
        deallog << ' ' << GeometryInfo<dim>::face_to_cell_vertices(f, v, true);
      deallog << std::endl;
      deallog << "face_vertices" << f << "[false]";
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
        deallog << ' ' << GeometryInfo<dim>::face_to_cell_vertices(f, v, false);
      deallog << std::endl;
    }

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    {
      deallog << "face_lines" << f << "[true ]";
      for (unsigned int v = 1; v <= GeometryInfo<dim>::lines_per_face; ++v)
        deallog << ' ' << GeometryInfo<dim>::face_to_cell_lines(f, v - 1, true);
      deallog << std::endl;
      deallog << "face_lines" << f << "[false]";
      for (unsigned int v = 1; v <= GeometryInfo<dim>::lines_per_face; ++v)
        deallog << ' '
                << GeometryInfo<dim>::face_to_cell_lines(f, v - 1, false);
      deallog << std::endl;
    }

  for (unsigned int f = 0; f < GeometryInfo<dim>::lines_per_cell; ++f)
    {
      deallog << "line_vertices" << f;
      for (const unsigned int v : GeometryInfo<1>::vertex_indices())
        deallog << ' ' << GeometryInfo<dim>::line_to_cell_vertices(f, v);
      deallog << std::endl;
    }
}


int
main()
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();
  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
  return 0;
}
