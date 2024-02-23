// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that face orientation flags work. at one point in time we
// mixed up the generation of edges associated with faces for the
// non-standard case when refining the mesh. like in mesh_3d_10, but
// make sure that each cell has 8 distinct vertices if looking at the
// vertices of the faces

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <set>

#include "../tests.h"

#include "mesh_3d.h"

void
check_this(Triangulation<3> &tria)
{
  for (Triangulation<3>::cell_iterator cell = tria.begin(); cell != tria.end();
       ++cell)
    {
      std::set<unsigned int> vertices;
      for (const unsigned int f : GeometryInfo<3>::face_indices())
        {
          vertices.insert(cell->face(f)->vertex_index(0));
          vertices.insert(cell->face(f)->vertex_index(1));
          vertices.insert(cell->face(f)->vertex_index(2));
          vertices.insert(cell->face(f)->vertex_index(3));
        }
      AssertThrow(vertices.size() == 8, ExcInternalError());
    }
  deallog << "    ok." << std::endl;
}


void
check(Triangulation<3> &tria)
{
  deallog << "Initial check" << std::endl;
  check_this(tria);

  for (unsigned int r = 0; r < 3; ++r)
    {
      tria.refine_global(1);
      deallog << "Check " << r << std::endl;
      check_this(tria);
    }

  coarsen_global(tria);
  deallog << "Check " << 1 << std::endl;
  check_this(tria);

  tria.refine_global(1);
  deallog << "Check " << 2 << std::endl;
  check_this(tria);
}


int
main()
{
  initlog();

  {
    Triangulation<3> coarse_grid;
    create_two_cubes(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_L_shape(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball(coarse_grid);
    check(coarse_grid);
  }
}
