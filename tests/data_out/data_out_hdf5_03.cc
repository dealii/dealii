// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test parallel DataOut with HDF5 + xdmf
//
// When running with 3 MPI ranks, one of the ranks will have 0
// cells. This tests a corner case that used to fail because the code
// assumed that all ranks have at least one cell.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include <string>
#include <vector>

#include "../tests.h"

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  FE_Q<dim> fe(1);

  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  Vector<double> v1(dof.n_dofs());
  for (unsigned int i = 0; i < v1.size(); ++i)
    v1(i) = i;

  DataOut<dim> data_out;
  data_out.add_data_vector(dof, v1, "bla");
  data_out.build_patches(1);

  DataOutBase::DataOutFilter data_filter(
    DataOutBase::DataOutFilterFlags(false, false));
  data_out.write_filtered_data(data_filter);
  deallog << "n_cells on my rank: " << data_filter.n_cells() << std::endl;
  data_out.write_hdf5_parallel(data_filter, "out.h5", MPI_COMM_WORLD);

  std::vector<XDMFEntry> xdmf_entries;
  xdmf_entries.push_back(
    data_out.create_xdmf_entry(data_filter, "out.h5", 0, MPI_COMM_WORLD));

  data_out.write_xdmf_file(xdmf_entries, "out.xdmf", MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      cat_file("out.xdmf");

      // Sadly hdf5 is binary and we can not use hd5dump because it might
      // not be in the path.
      std::ifstream f("out.h5");
      AssertThrow(f.good(), ExcIO());
      deallog << "ok" << std::endl;
    }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;


  test<2>();
  test<3>();
}
