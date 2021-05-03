// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2019 by the deal.II authors
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


// check the new DoFRenumbering::component_wise function that handles
// DoFHandlers and renumbers all MG and non-MG dofs

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <algorithm>

#include "../tests.h"

using namespace std;

template <int dim>
void
check()
{
  FESystem<dim> fe(FE_DGP<dim>(1), 1, FE_DGQ<dim>(2), 2, FE_Q<dim>(3), 1);
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  DoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  mgdof.distribute_mg_dofs();
  DoFRenumbering::component_wise(mgdof);
  for (unsigned int l = 0; l < tr.n_levels(); ++l)
    DoFRenumbering::component_wise(mgdof, l);


  typename DoFHandler<dim>::level_cell_iterator cell = mgdof.begin_mg(),
                                                endc = mgdof.end_mg();
  std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
  std::vector<types::global_dof_index> mg_dof_indices(fe.dofs_per_cell);
  for (; cell != endc; ++cell)
    {
      if (!cell->has_children())
        {
          deallog << "Global numbering: ";
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            deallog << local_dof_indices[i] << ' ';
          deallog << std::endl;
        }

      deallog << "MG levelwise numbering: ";
      cell->get_mg_dof_indices(mg_dof_indices);
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        deallog << mg_dof_indices[i] << ' ';
      deallog << std::endl;

      // assert at locally on each
      // cell that dofs with lower
      // component also have lower
      // dof index
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        for (unsigned int j = i + 1; j < fe.dofs_per_cell; ++j)
          {
            if (fe.system_to_component_index(i).first <
                fe.system_to_component_index(j).first)
              {
                AssertThrow(mg_dof_indices[i] < mg_dof_indices[j],
                            ExcInternalError());
              }
            else if (fe.system_to_component_index(i).first >
                     fe.system_to_component_index(j).first)
              {
                AssertThrow(mg_dof_indices[i] > mg_dof_indices[j],
                            ExcInternalError());
              }
          }
    }
}


int
main()
{
  initlog(__FILE__);
  check<1>();
  check<2>();
  check<3>();
}
