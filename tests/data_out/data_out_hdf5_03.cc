// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
//
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

// test parallel DataOut with HDF5
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
  std::vector<unsigned int>                 repetitions(dim, 1);
  repetitions[0] = 2; // 2x1x1 cells
  Point<dim> p1;
  Point<dim> p2;
  for (int i = 0; i < dim; ++i)
    p2[i] = 1.0;
  GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);
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
  data_out.write_hdf5_parallel(data_filter, "out.h5", MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
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
