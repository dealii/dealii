// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check GriOut::write_vtk() with subdomain ids

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"


template <int dim>
void
output(parallel::distributed::Triangulation<dim> &tr, std::ostream &stream)
{
  GridOut out;
  out.write_vtk(tr, stream);
}

template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  DoFHandler<dim> dofh(tr);

  // const std::string
  //   filename = ("mesh" + Utilities::int_to_string(dim)
  //    + "d." +
  //    Utilities::int_to_string(tr.locally_owned_subdomain(), 4) +
  //    ".vtk");
  // std::ofstream stream(filename);
  output(tr, deallog.get_file_stream());

  // tr.begin_active()->set_level_subdomain_id(1+myid);

  for (unsigned int i = 0; i < 10; ++i)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      if (myid == i)
        {
          deallog << "ID = " << i << std::endl;
          for (unsigned int lvl = 0; lvl < tr.n_levels(); ++lvl)
            {
              typename Triangulation<dim>::cell_iterator cell = tr.begin(lvl),
                                                         endc = tr.end(lvl);

              for (; cell != endc; ++cell)
                if (cell->level_subdomain_id() != 4294967294)
                  deallog << cell->level_subdomain_id();
                else
                  deallog << '-';
              deallog << std::endl;
            }
          deallog << std::endl;
        }
    }

  if (myid == 0)
    deallog << "my levels: " << tr.n_levels()
            << "<= global levels:" << tr.n_global_levels() << std::endl;
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
