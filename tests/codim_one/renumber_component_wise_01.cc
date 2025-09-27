// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test DoFRenumbering::component_wise for codim=1

#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

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

  FESystem<dim, spacedim>   fe(FE_Q<dim, spacedim>(2), spacedim);
  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;
  DoFRenumbering::component_wise(dof_handler);

  for (typename DoFHandler<dim, spacedim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      std::vector<types::global_dof_index> x(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(x);

      deallog << cell << std::endl;
      for (unsigned int i = 0; i < x.size(); ++i)
        deallog << "  " << x[i] << std::endl;
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
