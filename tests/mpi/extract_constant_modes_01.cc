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



// test DoFTools::extract_constant_modes

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
#include <deal.II/fe/fe_dgq.h>
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

  FESystem<dim> fe(FE_Q<dim>(1),2,
                   FE_DGQ<dim>(0),1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "Total dofs=" << dofh.n_dofs() << std::endl;

  // extract constant modes and print
  if (myid == 0)
    for (unsigned int c=0; c<fe.n_components(); ++c)
      {
        std::vector<bool> mask(fe.n_components(), false);
        mask[c] = true;

        std::vector<std::vector<bool> > constant_modes;
        DoFTools::extract_constant_modes (dofh, mask, constant_modes);

        for (unsigned int i=0; i<constant_modes.size(); ++i)
          {
            for (unsigned int j=0; j<constant_modes[i].size(); ++j)
              deallog << (constant_modes[i][j] ? '1' : '0') << ' ';
            deallog << std::endl;
          }
      }

  // renumber dofs and do the same again
  DoFRenumbering::component_wise(dofh);
  if (myid == 0)
    for (unsigned int c=0; c<fe.n_components(); ++c)
      {
        std::vector<bool> mask(fe.n_components(), false);
        mask[c] = true;

        std::vector<std::vector<bool> > constant_modes;
        DoFTools::extract_constant_modes (dofh, mask, constant_modes);

        for (unsigned int i=0; i<constant_modes.size(); ++i)
          {
            for (unsigned int j=0; j<constant_modes[i].size(); ++j)
              deallog << (constant_modes[i][j] ? '1' : '0') << ' ';
            deallog << std::endl;
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
