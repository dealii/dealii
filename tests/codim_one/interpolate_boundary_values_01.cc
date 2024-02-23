// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test VectorTools::interpolate_boundary_values for codim=1

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include <string>

#include "../tests.h"

template <int dim, int spacedim>
void
test(std::string filename)
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename);
  gi.read_ucd(in);

  deallog << tria.n_active_cells() << " active cells" << std::endl;

  FE_Q<dim, spacedim>       fe(2);
  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;

  std::map<types::global_dof_index, double> bv;
  VectorTools::interpolate_boundary_values(
    dof_handler, 0, Functions::SquareFunction<spacedim>(), bv);
  deallog << bv.size() << " boundary degrees of freedom" << std::endl;

  for (std::map<types::global_dof_index, double>::const_iterator i = bv.begin();
       i != bv.end();
       ++i)
    deallog << i->first << ' ' << i->second << std::endl;

  for (typename DoFHandler<dim, spacedim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
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



int
main()
{
  initlog();

  test<2, 3>(SOURCE_DIR "/grids/square.inp");
  test<2, 3>(SOURCE_DIR "/grids/sphere_1.inp");

  return 0;
}
