// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// Construct a sparsity pattern using NonMatching::DoFHandlerCoupling

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/non_matching/coupling.h>
#include <deal.II/non_matching/dof_handler_coupling.h>

#include <deal.II/numerics/matrix_tools.h>

#include "../tests.h"

using namespace dealii;

// Test that a coupling matrix can be constructed for each pair of dimension and
// immersed dimension, and check that constants are projected correctly.

template <int dim, int spacedim>
void
test()
{
  deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl;

  parallel::distributed::Triangulation<spacedim>      tria1(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim, spacedim> tria2(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria1, -1, 1);
  GridGenerator::hyper_cube(tria2, -.4, .3);

  tria1.refine_global(2);
  tria2.refine_global(1);

  FE_Q<spacedim>      fe1(1);
  FE_Q<dim, spacedim> fe2(1);

  deallog << "FE 1 : " << fe1.get_name() << std::endl
          << "FE 2 : " << fe2.get_name() << std::endl;

  DoFHandler<spacedim>      dh1(tria1);
  DoFHandler<dim, spacedim> dh2(tria2);

  dh1.distribute_dofs(fe1);
  dh2.distribute_dofs(fe2);

  deallog << "Dofs 1 : " << dh1.n_dofs() << std::endl
          << "Dofs 2 : " << dh2.n_dofs() << std::endl;

  NonMatching::DoFHandlerCoupling<dim, spacedim> dhc(dh1, dh2);

  const auto locally_owned2    = dh2.locally_owned_dofs();
  IndexSet   locally_relevant2 = DoFTools::extract_locally_relevant_dofs(dh2);
  const auto globally_owned2 =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, locally_owned2.n_elements());

  SparsityPattern sparsity;
  {
    DynamicSparsityPattern dsp(locally_relevant2);
    dhc.create_interpolation_sparsity_pattern(dsp);

    SparsityTools::distribute_sparsity_pattern(dsp,
                                               globally_owned2,
                                               MPI_COMM_WORLD,
                                               locally_relevant2);
    sparsity.copy_from(dsp);
  }
  sparsity.print(deallog.get_file_stream());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2, 2>();
}
