// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test Triangulation::communicate_locally_moved_vertices()

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  // create a mesh twice refined and move those vertices we locally own
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  const std::vector<bool> locally_owned_vertices =
    GridTools::get_locally_owned_vertices(tr);

  if (myid == 0)
    deallog << "#vertices = " << tr.n_vertices() << std::endl
            << "#locally_owned_vertices = "
            << std::count(locally_owned_vertices.begin(),
                          locally_owned_vertices.end(),
                          true)
            << std::endl;

  // now do the move
  Point<dim> shift;
  for (unsigned int d = 0; d < dim; ++d)
    shift[d] = 1;

  unsigned int n_vertices_moved = 0;
  for (unsigned int v = 0; v < tr.n_vertices(); ++v)
    if (locally_owned_vertices[v] == true)
      {
        // maybe not the most elegant way to do it, but it works for the purpose
        // of the test...
        const_cast<Point<dim> &>(tr.get_vertices()[v]) += shift;
        ++n_vertices_moved;
      }
  Assert(Utilities::MPI::sum(n_vertices_moved, MPI_COMM_WORLD) ==
           (dim == 2 ? 25 : 125),
         ExcInternalError());

  tr.communicate_locally_moved_vertices(locally_owned_vertices);

  // let every processor write their triangulation into a file of their own,
  // then collate
  // write the info on ghost processors and import indices to file
  {
    std::ofstream file((std::string("communicate_moved_vertices_01.dat.") +
                        Utilities::int_to_string(myid))
                         .c_str());
    GridOut().write_gnuplot(tr, file);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid == 0)
    {
      for (unsigned int i = 0;
           i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
           ++i)
        {
          deallog << "Partition " << i << std::endl;

          cat_file((std::string("communicate_moved_vertices_01.dat.") +
                    Utilities::int_to_string(i))
                     .c_str());
        }
    }


  if (myid == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
