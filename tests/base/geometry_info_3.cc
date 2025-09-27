// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check GeometryInfo::face_to_cell_vertices

#include <deal.II/base/geometry_info.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "Checking in " << dim << 'd' << std::endl;

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
