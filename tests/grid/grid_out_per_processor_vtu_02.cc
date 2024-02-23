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


// check GriOut::write_mesh_per_processor_as_vtu() when level_subdomain_id
// differs from subdomain_id for a particular cell

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>
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

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
void
output(const parallel::shared::Triangulation<dim> &tr,
       const std::string                          &filename,
       const bool                                  view_levels,
       const bool                                  include_artificial)
{
  GridOut           out;
  GridOutFlags::Vtu vtu_flags;
  vtu_flags.compression_level = DataOutBase::CompressionLevel::best_compression;
  out.set_flags(vtu_flags);
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

  parallel::shared::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_metis);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  typename Triangulation<dim>::active_cell_iterator cell = tr.begin_active(),
                                                    endc = tr.end();
  for (; cell != endc; ++cell)
    {
      if (cell->index() < 2)
        cell->set_subdomain_id(cell->index());
      else
        cell->set_subdomain_id(numbers::artificial_subdomain_id);

      if (cell->index() == 0 || cell->index() == 2)
        cell->set_level_subdomain_id(numbers::artificial_subdomain_id);
      else if (cell->index() == 1)
        cell->set_level_subdomain_id(0);
      else if (cell->index() == 3)
        cell->set_level_subdomain_id(1);
    }

  output(tr, "file1", true, false);
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
