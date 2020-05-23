// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
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


// check DoFAccessor::get_mg_dof_indices on meshes that are not in standard
// orientation by comparing to DoFAccessor::get_dof_indices

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
check()
{
  // need cubic polynomials that have two dofs on lines
  FE_Q<dim> fe(3);

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  if (dim > 1)
    GridGenerator::hyper_shell(tr, Point<dim>(), 0.5, 1, 12);
  else
    GridGenerator::hyper_cube(tr);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs();

  std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
  std::vector<types::global_dof_index> mg_dof_indices(fe.dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end();
       ++cell)
    {
      cell->get_dof_indices(dof_indices);
      cell->get_mg_dof_indices(mg_dof_indices);
      bool has_error = false;
      // dof indices should have the same order on both the mg dofs and the
      // usual dofs because there is only one level
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        if (dof_indices[i] != mg_dof_indices[i])
          has_error = true;
      if (has_error)
        {
          deallog << "Offending cell with center " << cell->center() << ": ";
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            deallog << dof_indices[i] << " vs " << mg_dof_indices[i] << ", ";
          deallog << std::endl;
          return;
        }
    }
  deallog << dim << "D OK" << std::endl;
}

int
main()
{
  initlog(__FILE__);
  check<1>();
  check<2>();
  check<3>();
}
