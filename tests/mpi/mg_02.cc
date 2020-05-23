// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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



// distribute dofs and mgdofs for a parallel Triangulation

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  DoFHandler<dim> dofh(tr);

  {
    for (unsigned int lvl = 0; lvl < tr.n_levels(); ++lvl)
      {
        deallog << "level " << lvl << ": ";
        typename DoFHandler<dim>::cell_iterator cell = dofh.begin(lvl),
                                                endc = dofh.end(lvl);

        for (; cell != endc; ++cell)
          {
            if (cell->level_subdomain_id() != 4294967294)
              deallog << cell->level_subdomain_id();
            else
              deallog << "-";
          }
        deallog << std::endl;
      }
  }

  static const FE_DGP<dim> fe(0);
  Assert(dofh.has_active_dofs() == false, ExcInternalError());
  Assert(dofh.has_level_dofs() == false, ExcInternalError());

  dofh.distribute_dofs(fe);

  Assert(dofh.has_active_dofs() == true, ExcInternalError());
  Assert(dofh.has_level_dofs() == false, ExcInternalError());

  dofh.distribute_mg_dofs();

  Assert(dofh.has_active_dofs() == true, ExcInternalError());
  Assert(dofh.has_level_dofs() == true, ExcInternalError());
  {
    deallog << "Levels: " << tr.n_global_levels() << std::endl;
    std::cout << "Levels: " << tr.n_global_levels() << std::endl;

    const std::vector<types::global_dof_index>
      n_locally_owned_dofs_per_processor =
        Utilities::MPI::all_gather(MPI_COMM_WORLD, dofh.n_locally_owned_dofs());
    deallog << "n_locally_owned_dofs_per_processor:" << std::endl;
    for (unsigned int i = 0; i < n_locally_owned_dofs_per_processor.size(); ++i)
      deallog << n_locally_owned_dofs_per_processor[i] << std::endl;

    deallog << "locally_owned_mg_dofs_per_processor:" << std::endl;
    for (unsigned int lvl = 0; lvl < tr.n_global_levels(); ++lvl)
      {
        deallog << "level " << lvl << ":" << std::endl;

        const std::vector<IndexSet> vec =
          Utilities::MPI::all_gather(MPI_COMM_WORLD,
                                     dofh.locally_owned_mg_dofs(lvl));

        for (unsigned int i = 0; i < vec.size(); ++i)
          deallog << vec[i].n_elements() << std::endl;
      }
  }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();
}
