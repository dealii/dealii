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



// create a parallel DoFHandler on a 2d mesh and check componentwise
// renumbering
//
// this is like the _01 testcase but with a simpler mesh for which I
// can actually count DoFs and decide who owns what


#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>


template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  std::vector<unsigned int> sub(2);
  sub[0] = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);
  sub[1] = 1;
  GridGenerator::subdivided_hyper_rectangle(static_cast<Triangulation<dim>&>(tr),
                                            sub, Point<2>(0,0), Point<2>(1,1));

  const FE_Q<dim> fe_q(1);
  FESystem<dim> fe(fe_q,2);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "Total dofs=" << dofh.n_dofs() << std::endl;

  {
    IndexSet dof_set;
    DoFTools::extract_locally_active_dofs (dofh, dof_set);
    if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
      dof_set.print(deallog);
  }

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "****" << std::endl;

  DoFRenumbering::component_wise(dofh);
  {
    IndexSet dof_set;
    DoFTools::extract_locally_active_dofs (dofh, dof_set);

    if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
      {
        dof_set.print(deallog);
        for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
          {
            deallog << "Dofs owned by processor " << i << ": ";
            dofh.locally_owned_dofs_per_processor()[i].print(deallog);
            deallog << std::endl;
          }
      }

    if (myid==0)
      {

        std::vector<types::global_dof_index> local_dof_indices;
        typename DoFHandler<dim>::active_cell_iterator
        cell, endc = dofh.end();

        if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
          for (cell = dofh.begin_active(); cell != endc; ++cell)
            if (!cell->is_artificial() && !cell->is_ghost())
              {
                local_dof_indices.resize (cell->get_fe().dofs_per_cell);
                cell->get_dof_indices (local_dof_indices);
                for (unsigned int i=0; i< cell->get_fe().dofs_per_cell; ++i)
                  deallog << local_dof_indices[i] << " ";
                deallog << std::endl;
              }
      }
  }
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();
    }

}
