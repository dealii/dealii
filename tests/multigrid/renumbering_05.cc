// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// check renumbering the degrees of freedom on the multigrid levels in
// parallel

#include "../tests.h"
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <algorithm>

using namespace std;

template <int dim>
void
print_dof_numbers(const DoFHandler<dim>& dof)
{
  std::vector<types::global_dof_index> dof_indices(dof.get_fe().dofs_per_cell);
  deallog << "DoF numbers on active cells" << std::endl;
  for(auto cell : dof.active_cell_iterators())
    if(!cell->is_artificial())
      {
        cell->get_dof_indices(dof_indices);
        deallog << "cell " << cell->id() << ": ";
        for(types::global_dof_index i : dof_indices)
          deallog << i << " ";
        deallog << std::endl;
      }
  for(unsigned int l = 0; l < dof.get_triangulation().n_global_levels(); ++l)
    {
      deallog << "DoF numbers on level " << l << std::endl;
      for(auto cell : dof.cell_iterators_on_level(l))
        if(cell->level_subdomain_id() != numbers::artificial_subdomain_id)
          {
            cell->get_mg_dof_indices(dof_indices);
            deallog << "cell " << cell->id() << ": ";
            for(types::global_dof_index i : dof_indices)
              deallog << i << " ";
            deallog << std::endl;
          }
    }
}

template <int dim>
void
check()
{
  FE_DGQ<dim> fe(1);

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(5 - dim);
  if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    tr.begin_active(tr.n_global_levels() - 1)->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  DoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  mgdof.distribute_mg_dofs();

  print_dof_numbers(mgdof);

  // compute a renumbering on the active degrees of freedom
  {
    std::vector<types::global_dof_index> new_indices;
    if(mgdof.n_locally_owned_dofs() > 0)
      {
        const types::global_dof_index first
          = mgdof.locally_owned_dofs().nth_index_in_set(0);
        const types::global_dof_index last
          = first + mgdof.n_locally_owned_dofs();
        const unsigned int stride = (last - first) / fe.dofs_per_cell;
        for(unsigned int j = 0; j < stride; ++j)
          for(unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            new_indices.push_back(first + j + i * stride);
      }
    mgdof.renumber_dofs(new_indices);
  }

  // compute a renumbering on the level degrees of freedom
  for(unsigned int l = 0; l < tr.n_global_levels(); ++l)
    {
      std::vector<types::global_dof_index> new_indices;
      if(mgdof.locally_owned_mg_dofs(l).n_elements() > 0)
        {
          const types::global_dof_index first
            = mgdof.locally_owned_mg_dofs(l).nth_index_in_set(0);
          const types::global_dof_index last
            = first + mgdof.locally_owned_mg_dofs(l).n_elements();
          const unsigned int stride = (last - first) / fe.dofs_per_cell;
          for(unsigned int j = 0; j < stride; ++j)
            for(unsigned int i = 0; i < fe.dofs_per_cell; ++i)
              new_indices.push_back(first + j + i * stride);
        }
      mgdof.renumber_dofs(l, new_indices);
    }

  print_dof_numbers(mgdof);
}

int
main(int argc, char** argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);
  MPILogInitAll                    log;

  check<2>();
  check<3>();
}
