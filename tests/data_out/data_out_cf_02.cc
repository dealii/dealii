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

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

#include "./data_out_cf_common.h"

template <int dim>
void
test()
{
  MPI_Comm                                  comm = MPI_COMM_WORLD;
  parallel::distributed::Triangulation<dim> tria(comm);
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);

  FE_Q<dim> fe1(1);

  DoFHandler<dim> dof1(tria);
  dof1.distribute_dofs(fe1);

  for (unsigned int t_index = 0; t_index < 3; ++t_index)
    {
      Vector<double> v1(dof1.n_dofs());
      for (unsigned int i = 0; i < v1.size(); ++i)
        v1(i) = i + t_index;

      DataOutBase::CFFlags flags;
      flags.time               = t_index * 1.0;
      flags.keep_existing_file = (t_index > 0);
      flags.attributes["global"].push_back({"comment", "test data_out_cf_02"});
      flags.attributes["time"].push_back({"units", "s since 1970-1-1"});
      flags.attributes["linear"].push_back(
        {"comment", "Interpolated DoF index."});

      DataOut<dim> data_out;
      data_out.set_flags(flags);
      data_out.add_data_vector(dof1, v1, "linear");
      data_out.build_patches(/*subdivisions*/ 2);

      DataOutBase::DataOutFilter data_filter(
        DataOutBase::DataOutFilterFlags(/*filter_vertices*/ true));
      data_out.write_filtered_data(data_filter);
      data_out.write_cf_parallel(data_filter, "out.nc", comm);
    }

  if (Utilities::MPI::this_mpi_process(comm) == 0)
    {
      deallog << "ok" << std::endl;
      dump_nc_file(deallog, "out.nc");
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

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
