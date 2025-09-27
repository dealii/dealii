// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test VectorTools::interpolate_boundary_values for codim=1. like _01
// but for 1d triangulations

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include <string>

#include "../tests.h"

void
test()
{
  const int dim      = 1;
  const int spacedim = 2;

  Triangulation<dim, spacedim> tria;
  Triangulation<spacedim>      volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);
  const std::set<types::boundary_id> boundary_ids = {0};
  GridGenerator::extract_boundary_mesh(volume_mesh, tria, boundary_ids);

  deallog << tria.n_active_cells() << " active cells" << std::endl;

  FE_Q<dim, spacedim>       fe(1);
  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;

  // test left and right boundary
  // separately
  for (unsigned int boundary_id = 0; boundary_id < 2; ++boundary_id)
    {
      std::map<types::global_dof_index, double> bv;
      VectorTools::interpolate_boundary_values(
        dof_handler, boundary_id, Functions::SquareFunction<spacedim>(), bv);
      deallog << bv.size() << " boundary degrees of freedom" << std::endl;

      for (std::map<types::global_dof_index, double>::const_iterator i =
             bv.begin();
           i != bv.end();
           ++i)
        deallog << i->first << ' ' << i->second << std::endl;

      for (DoFHandler<dim, spacedim>::active_cell_iterator cell =
             dof_handler.begin_active();
           cell != dof_handler.end();
           ++cell)
        for (const unsigned int f : GeometryInfo<dim>::face_indices())
          if (cell->at_boundary(f) &&
              (cell->face(f)->boundary_id() == boundary_id))
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
                 ++v)
              for (unsigned int i = 0; i < fe.dofs_per_vertex; ++i)
                {
                  AssertThrow(bv.find(cell->face(f)->vertex_dof_index(v, i)) !=
                                bv.end(),
                              ExcInternalError());
                  AssertThrow(bv[cell->face(f)->vertex_dof_index(v, i)] ==
                                Functions::SquareFunction<spacedim>().value(
                                  cell->face(f)->vertex(v), i),
                              ExcInternalError());
                }
    }
}



int
main()
{
  initlog();

  test();

  return 0;
}
