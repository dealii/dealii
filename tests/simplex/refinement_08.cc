// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Distribute FE_WedgeP on a DoFHandler while refining the wedge (similar to
// wedge_01)

#include "deal.II/base/logstream.h"
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "./simplex_grids.h"


void
test()
{
  constexpr int      dim = 3;
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_wedges(tria, 1);

  // total number of refinements
  const unsigned int n_refinement = 2;

  for (unsigned int i = 0; i < n_refinement; i++)
    {
      FE_WedgeP<dim>  fe(1);
      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      deallog << "Number of refinements: " << i << std::endl;

      // print out the dof indices (see test wedge_01)
      std::vector<types::global_dof_index> dof_indices;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          for (const auto face : cell->face_indices())
            {
              dof_indices.resize(fe.n_dofs_per_face(face));
              cell->face(face)->get_dof_indices(dof_indices);

              for (const auto i : dof_indices)
                deallog << i << ' ';
              deallog << std::endl;
            }
        }

      // refine once for the next step
      if (i < n_refinement - 1)
        tria.refine_global(1);
    }
}

int
main()
{
  initlog();
  test();
}
