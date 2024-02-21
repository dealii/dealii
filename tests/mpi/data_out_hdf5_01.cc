// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// tests DataOut with HDF5

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);

  FE_Q<dim> fe1(1);

  DoFHandler<dim> dof1(tria);
  dof1.distribute_dofs(fe1);

  Vector<double> v1(dof1.n_dofs());
  for (unsigned int i = 0; i < v1.size(); ++i)
    v1(i) = i;

  DataOut<dim> data_out;
  data_out.add_data_vector(dof1, v1, "linear");
  data_out.build_patches(2);

  DataOutBase::DataOutFilter data_filter(
    DataOutBase::DataOutFilterFlags(false, false));
  data_out.write_filtered_data(data_filter);
  data_out.write_hdf5_parallel(data_filter, "out.h5", MPI_COMM_WORLD);
  std::vector<XDMFEntry> xdmf_entries;
  xdmf_entries.push_back(
    data_out.create_xdmf_entry(data_filter, "out.h5", 0, MPI_COMM_WORLD));

  data_out.write_xdmf_file(xdmf_entries, "out.xdmf", MPI_COMM_WORLD);

  deallog << "ok" << std::endl;


  // Sadly hdf5 is binary and we can not use hd5dump because it might
  // not be in the path. At least we can look at the xdmf
  // and make sure that the h5 file is created:
  if (0 == Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    {
      cat_file("out.xdmf");
      std::ifstream f("out.h5");
      AssertThrow(f.good(), ExcIO());
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  try
    {
      test<2>();

      return 0;
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
