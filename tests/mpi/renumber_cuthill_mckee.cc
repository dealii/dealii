// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// Test that DofRenumbering::Cuthill_McKee works in parallel (by applying
// Cuthill-McKee individually on each processor's subdomain).


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int nprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr (MPI_COMM_WORLD);

  GridGenerator::hyper_cube (tr, -1.0, 1.0);
  tr.refine_global (8-2*dim);

  for (typename Triangulation<dim>::active_cell_iterator
       cell = tr.begin_active();
       cell != tr.end(); ++cell)
    if (!cell->is_ghost() && !cell->is_artificial())
      if (cell->center().norm() < 0.3)
        {
          cell->set_refine_flag();
        }

  tr.execute_coarsening_and_refinement ();

  DoFHandler<dim> dofh (tr);

  static const FE_Q<dim> fe (1);
  dofh.distribute_dofs (fe);
  std::vector<types::global_dof_index> renumbering(dofh.locally_owned_dofs().n_elements());
  DoFRenumbering::compute_Cuthill_McKee (renumbering, dofh);

  // send everything to processor 0 for output
  std::vector<types::global_dof_index> complete_renumbering(dofh.n_dofs());
  std::copy(renumbering.begin(), renumbering.end(), complete_renumbering.begin());
  unsigned int offset = renumbering.size();
  for (unsigned int i=1; i<nprocs; ++i)
    {
      if (myid == i)
        MPI_Send (&renumbering[0], renumbering.size(),
                  Utilities::MPI::internal::mpi_type_id(&complete_renumbering[0]), 0, i,
                  MPI_COMM_WORLD);
      else if (myid == 0)
        MPI_Recv (&complete_renumbering[offset],
                  dofh.locally_owned_dofs_per_processor()[i].n_elements(),
                  Utilities::MPI::internal::mpi_type_id(&complete_renumbering[0]), i, i,
                  MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
      offset += dofh.locally_owned_dofs_per_processor()[i].n_elements();
    }

  if (myid == 0)
    {
      AssertDimension(offset, complete_renumbering.size());
      for (unsigned int i=0; i<complete_renumbering.size(); ++i)
        deallog << complete_renumbering[i] << std::endl;
    }
}


int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push (Utilities::int_to_string (myid));

  if (myid == 0)
    {
      std::ofstream logfile ("output");
      deallog.attach (logfile);
      deallog.depth_console (0);
      deallog.threshold_double (1.e-10);

      deallog.push ("2d");
      test<2>();
      deallog.pop();
      deallog.push ("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }
}
