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



// check GriOut::write_mesh_per_processor_as_vtu()

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
void
output(const Triangulation<dim> &tr,
       const std::string &       filename,
       const bool                view_levels,
       const bool                include_artificial)
{
  GridOut out;
  out.write_mesh_per_processor_as_vtu(tr,
                                      filename,
                                      view_levels,
                                      include_artificial);

  // copy the .pvtu and .vtu files
  // into the logstream
  int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (myid == 0)
    {
      cat_file((std::string(filename) + ".pvtu").c_str());
      cat_file((std::string(filename) + ".proc0000.vtu").c_str());
    }
  else if (myid == 1)
    cat_file((std::string(filename) + ".proc0001.vtu").c_str());
  else if (myid == 2)
    cat_file((std::string(filename) + ".proc0002.vtu").c_str());
  else
    AssertThrow(false, ExcNotImplemented());
}

template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "hyper_cube" << std::endl;


  parallel::shared::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::Settings::partition_zorder);

  Triangulation<dim> temp;
  GridGenerator::hyper_cube(temp);
  GridGenerator::convert_hypercube_to_simplex_mesh(temp, triangulation);

  triangulation.refine_global(2);

  DoFHandler<dim> dofh(triangulation);

  output(triangulation, "file1", false, true);
  output(triangulation, "file2", true, false);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll init;

  deallog.push("2d");
  test<2>();
  deallog.pop();
}
