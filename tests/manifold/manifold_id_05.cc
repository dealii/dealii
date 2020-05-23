// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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

// Set a manifold id on one of the boundary face, and attach the
// boundary description to it. Refine globally twice and output the mesh.

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// Helper function
template <int dim, int spacedim>
void
test(unsigned int ref = 1)
{
  deallog << "Testing dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Point<spacedim> center;
  for (unsigned int i = 0; i < dim; ++i)
    center[i] = .5;

  SphericalManifold<dim, spacedim> boundary(center);
  Triangulation<dim, spacedim>     tria;
  GridGenerator::hyper_cube(tria);
  typename Triangulation<dim, spacedim>::active_cell_iterator cell;

  tria.begin_active()->face(0)->set_manifold_id(1);
  tria.set_manifold(1, boundary);

  tria.refine_global(2);

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      deallog << "C: " << cell << std::endl;
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        deallog << "F: " << cell->face(f)
                << ", mid: " << (int)cell->face(f)->manifold_id() << std::endl;
    }


  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  return 0;
}
