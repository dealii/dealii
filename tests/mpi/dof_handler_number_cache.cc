// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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



// like the same test in deal.II but this time use a
// parallel::distributed::Triangulation object. We still use only a
// single processor so the end result should be the same but we use
// entirely different code paths


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>

#include <fstream>
#include <cstdlib>
#include <numeric>


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim>
  triangulation (MPI_COMM_WORLD,
                 Triangulation<dim>::limit_level_difference_at_vertices);

  FESystem<dim> fe (FE_Q<dim>(3),2,
                    FE_DGQ<dim>(1),1);

  DoFHandler<dim> dof_handler (triangulation);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global (2);

  const unsigned int n_refinements[] = { 0, 4, 3, 2 };
  for (unsigned int i=0; i<n_refinements[dim]; ++i)
    {
      // refine one-fifth of cells randomly
      std::vector<bool> flags (triangulation.n_active_cells(), false);
      for (unsigned int k=0; k<flags.size()/5 + 1; ++k)
        flags[Testing::rand() % flags.size()] = true;
      // make sure there's at least one that
      // will be refined
      flags[0] = true;

      // refine triangulation
      unsigned int index=0;
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell)
        if (!cell->is_ghost() && !cell->is_artificial())
          {
            if (flags[index])
              cell->set_refine_flag();
            ++index;
          }

      Assert (index <= triangulation.n_active_cells(), ExcInternalError());

      // flag all other cells for coarsening
      // (this should ensure that at least
      // some of them will actually be
      // coarsened)
      index=0;
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell)
        if (!cell->is_ghost() && !cell->is_artificial())
          {
            if (!flags[index])
              cell->set_coarsen_flag();
            ++index;
          }

      triangulation.execute_coarsening_and_refinement ();
      dof_handler.distribute_dofs (fe);

      const unsigned int N = dof_handler.n_dofs();
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        deallog << N << std::endl;

      Assert (dof_handler.n_locally_owned_dofs() <= N,
              ExcInternalError());
      for (unsigned int i=0;
           i<dof_handler.n_locally_owned_dofs_per_processor().size(); ++i)
        Assert (dof_handler.n_locally_owned_dofs_per_processor()[i] <= N,
                ExcInternalError());
      Assert (std::accumulate (dof_handler.n_locally_owned_dofs_per_processor().begin(),
                               dof_handler.n_locally_owned_dofs_per_processor().end(),
                               0U) == N,
              ExcInternalError());

      IndexSet all (N), really_all (N);
      // poor man's union operation
      for (unsigned int i=0;
           i<dof_handler.n_locally_owned_dofs_per_processor().size(); ++i)
        for (unsigned int j=0; j<N; ++j)
          if (dof_handler.locally_owned_dofs_per_processor()[i]
              .is_element(j))
            {
              Assert (all.is_element(j) == false, ExcInternalError());
              all.add_index (j);
            }
      really_all.add_range (0, N);
      Assert (all == really_all,
              ExcInternalError());
    }
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }

}
