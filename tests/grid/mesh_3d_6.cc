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



// check that face orientation flags work by looping over all cells
// and check on all faces that if we look from both sides that normal
// vectors point in opposite directions

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "mesh_3d.h"



void
check_this(Triangulation<3> &tria)
{
  QMidpoint<2>    q;
  FE_Q<3>         fe(1);
  FEFaceValues<3> fe_face_values1(fe, q, update_normal_vectors);
  FEFaceValues<3> fe_face_values2(fe, q, update_normal_vectors);

  DoFHandler<3> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  unsigned int global_face = 0;

  // look at all faces, not only
  // active ones
  for (DoFHandler<3>::cell_iterator cell = dof_handler.begin();
       cell != dof_handler.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<3>::face_indices())
      if (!cell->at_boundary(f))
        {
          const unsigned int nn = cell->neighbor_of_neighbor(f);
          fe_face_values1.reinit(cell, f);
          fe_face_values2.reinit(cell->neighbor(f), nn);

          // in order to reduce
          // output file size, only
          // write every seventeenth
          // normal vector. if the
          // normals differ anyway,
          // then the assertion below
          // will catch this, and if
          // we compute _all_ normals
          // wrongly, then outputting
          // some will be ok, I guess
          if (global_face++ % 17 == 0)
            deallog << "Cell " << cell << ", face " << f
                    << " n=" << fe_face_values1.normal_vector(0) << std::endl;

          // normal vectors should be
          // in opposite directions,
          // so their sum should be
          // close to zero
          Assert((fe_face_values1.normal_vector(0) +
                  fe_face_values2.normal_vector(0))
                     .norm_square() < 1e-20,
                 ExcInternalError());
        }
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
    coarse_grid.reset_manifold(0);
    check(coarse_grid);
  }
}
