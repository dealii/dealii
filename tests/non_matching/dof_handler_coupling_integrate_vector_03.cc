// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Integrate a constant field represented on dh1 against basis functions on dh2.

#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/non_matching/dof_handler_coupling.h>

#include "../tests.h"

using namespace dealii;


template <int dim, int spacedim>
void
test()
{
  deallog << "integrate_vector_01" << std::endl;

  parallel::distributed::Triangulation<spacedim>      tria1(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim, spacedim> tria2(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria1, -1.0, 1.0);
  GridGenerator::hyper_cube(tria2, -0.8, 0.7);

  tria1.refine_global(3);
  tria2.refine_global(2);

  FE_Q<spacedim>      fe1(1);
  FE_Q<dim, spacedim> fe2(1);

  DoFHandler<spacedim>      dh1(tria1);
  DoFHandler<dim, spacedim> dh2(tria2);

  dh1.distribute_dofs(fe1);
  dh2.distribute_dofs(fe2);

  deallog << "n_dofs(dh1): " << dh1.n_dofs() << std::endl;
  deallog << "n_dofs(dh2): " << dh2.n_dofs() << std::endl;

  NonMatching::DoFHandlerCoupling<dim, spacedim> dhc(dh1, dh2);

  const IndexSet locally_owned_1 = dh1.locally_owned_dofs();
  const IndexSet locally_owned_2 = dh2.locally_owned_dofs();
  IndexSet locally_relevant_1    = DoFTools::extract_locally_relevant_dofs(dh1);

  LinearAlgebra::distributed::Vector<double> src;
  src.reinit(locally_owned_1, locally_relevant_1, MPI_COMM_WORLD);
  src=1.;
  src.update_ghost_values();

  QGauss<dim> quadrature(2);

  IndexSet locally_relevant_2 = dhc.extract_immersed_dof_indexset(quadrature);

  LinearAlgebra::distributed::Vector<double> dst;
  dst.reinit(locally_owned_2, locally_relevant_2, MPI_COMM_WORLD);

  dhc.integrate_dh1_field_against_dh2_basis(quadrature, src, dst);

  const double l1 = dst.l1_norm();
  const double l2 = dst.l2_norm();

  deallog << "dst l1_norm: " << l1 << std::endl;
  deallog << "dst l2_norm: " << l2 << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test<2, 2>();
}
