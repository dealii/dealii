// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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



// check number cache for shared_tria

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <numeric>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_metis);

  FESystem<dim> fe(FE_Q<dim>(3), 2, FE_DGQ<dim>(1), 1);

  DoFHandler<dim> dof_handler(triangulation);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  const unsigned int n_refinements[] = {0, 4, 3, 2};
  for (unsigned int i = 0; i < n_refinements[dim]; ++i)
    {
      // refine one-fifth of cells randomly
      std::vector<bool> flags(triangulation.n_active_cells(), false);
      for (unsigned int k = 0; k < flags.size() / 5 + 1; ++k)
        flags[Testing::rand() % flags.size()] = true;
      // make sure there's at least one that
      // will be refined
      flags[0] = true;

      // refine triangulation
      unsigned int index = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        {
          if (flags[index])
            cell->set_refine_flag();
          ++index;
        }

      Assert(index <= triangulation.n_active_cells(), ExcInternalError());

      // flag all other cells for coarsening
      // (this should ensure that at least
      // some of them will actually be
      // coarsened)
      index = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        {
          if (!flags[index])
            cell->set_coarsen_flag();
          ++index;
        }

      triangulation.execute_coarsening_and_refinement();
      dof_handler.distribute_dofs(fe);

      // avoid outputting any partitioning info until Parmetis gives
      // reproducible results
      deallog << "n_dofs: " << dof_handler.n_dofs() << std::endl;
      //          << "n_locally_owned_dofs: " <<
      //          dof_handler.n_locally_owned_dofs() << std::endl;

      //      deallog << "n_locally_owned_dofs_per_processor: ";
      //      std::vector<types::global_dof_index> v =
      //        dof_handler.compute_n_locally_owned_dofs_per_processor();
      //      unsigned int sum = 0;
      //      for (unsigned int i=0; i<v.size(); ++i)
      //        {
      //          deallog << v[i] << " ";
      //          sum += v[i];
      //        }
      //      deallog << " sum: " << sum << std::endl;

      const std::vector<types::global_dof_index>
        n_locally_owned_dofs_per_processor =
          Utilities::MPI::all_gather(MPI_COMM_WORLD,
                                     dof_handler.n_locally_owned_dofs());
      Assert(dof_handler.n_locally_owned_dofs() ==
               n_locally_owned_dofs_per_processor[triangulation
                                                    .locally_owned_subdomain()],
             ExcInternalError());
      Assert(dof_handler.n_locally_owned_dofs() ==
               dof_handler.locally_owned_dofs().n_elements(),
             ExcInternalError());

      const unsigned int N = dof_handler.n_dofs();

      Assert(dof_handler.n_locally_owned_dofs() <= N, ExcInternalError());
      Assert(std::accumulate(n_locally_owned_dofs_per_processor.begin(),
                             n_locally_owned_dofs_per_processor.end(),
                             0U) == N,
             ExcInternalError());

      const std::vector<IndexSet> locally_owned_dofs_per_processor =
        Utilities::MPI::all_gather(MPI_COMM_WORLD,
                                   dof_handler.locally_owned_dofs());
      IndexSet all(N);
      for (unsigned int i = 0; i < locally_owned_dofs_per_processor.size(); ++i)
        {
          IndexSet intersect = all & locally_owned_dofs_per_processor[i];
          Assert(intersect.n_elements() == 0, ExcInternalError());
          all.add_indices(locally_owned_dofs_per_processor[i]);
        }

      Assert(all == complete_index_set(N), ExcInternalError());
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
