// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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



// check GeometryInfo::face_to_cell_vertices

#include <deal.II/base/geometry_info.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "Checking in " << dim << "d" << std::endl;

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
      {
        deallog << "Face " << f << ", vertex=" << v << ": ";
        deallog << GeometryInfo<dim>::face_to_cell_vertices(f, v, true)
                << std::endl;
      }

  if (dim == 3)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
        {
          deallog << "Face " << f << ", vertex=" << v
                  << " (reverse orientation): ";
          deallog << GeometryInfo<dim>::face_to_cell_vertices(f, v, false)
                  << std::endl;
        }
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
