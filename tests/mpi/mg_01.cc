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



// check level_subdomain_id for distributed Triangulation

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/trilinos_vector.h>

#include <fstream>

template<int dim>
void output(parallel::distributed::Triangulation<dim> &tr)
{
  const std::string filename = ("mesh." +
                                Utilities::int_to_string
                                (tr.locally_owned_subdomain(), 4) +
                                ".fig");
  std::ofstream stream(filename.c_str());

  GridOutFlags::XFig flags;
  flags.color_by = GridOutFlags::XFig::level_subdomain_id;
  GridOut out;
  out.set_flags(flags);

  out.write_xfig(tr, stream);
}

template<int dim>
void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
                                               Triangulation<dim>::none,
                                               parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  DoFHandler<dim> dofh(tr);

  output(tr);

  //tr.begin_active()->set_level_subdomain_id(1+myid);

  for (unsigned int i=0; i<10; ++i)
    {

      MPI_Barrier(MPI_COMM_WORLD);
      if (myid==i)
        {
          deallog << "ID = " << i << std::endl;
          for (unsigned int lvl=0; lvl<tr.n_levels(); ++lvl)
            {
              typename Triangulation<dim>::cell_iterator
              cell = tr.begin(lvl),
              endc = tr.end(lvl);

              for (; cell!=endc; ++cell)
                if (cell->level_subdomain_id()!=4294967294)
                  deallog << cell->level_subdomain_id();
                else
                  deallog << "-";
              deallog << std::endl;
            }
          deallog << std::endl;
        }
    }

  if (myid==0)
    deallog << "my levels: " << tr.n_levels() << "<= global levels:" << tr.n_global_levels() << std::endl;


}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  deallog.push("2d");
  test<2>();
  deallog.pop();
}
