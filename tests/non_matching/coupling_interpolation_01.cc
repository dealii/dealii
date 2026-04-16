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

// Test that a coupling mass matrix can be constructed using
// DoFHandlerCoupling for each pair of dimension and immersed dimension,
// and check that constants are projected correctly.
// This is equivalent to coupling_01.cc but uses DoFHandlerCoupling instead of
// the free functions in non_matching/coupling.h.

#include <deal.II/base/mpi_stub.h>
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

#include <deal.II/non_matching/dof_handler_coupling.h>

#include <deal.II/numerics/matrix_tools.h>

#include "../tests.h"


// Test that a coupling mass matrix can be constructed for each pair of dimension
// and immersed dimension, and check that constants are projected correctly.

template <int dim, int spacedim>
void
test()
{
  deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl;

  parallel::distributed::Triangulation<dim, spacedim> tria(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<spacedim>      space_tria(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria, -.4, .3);
  GridGenerator::hyper_cube(space_tria, -1, 1);

  tria.refine_global(1);
  space_tria.refine_global(2);

  FE_Q<dim, spacedim> fe(1);
  FE_Q<spacedim>      space_fe(1);

  deallog << "FE      : " << fe.get_name() << std::endl
          << "Space FE: " << space_fe.get_name() << std::endl;

  DoFHandler<dim, spacedim> dh(tria);
  DoFHandler<spacedim>      space_dh(space_tria);

  dh.distribute_dofs(fe);
  space_dh.distribute_dofs(space_fe);

  deallog << "Dofs      : " << dh.n_dofs() << std::endl
          << "Space dofs: " << space_dh.n_dofs() << std::endl;

  QGauss<dim> quad(3); // Quadrature for coupling

  NonMatching::DoFHandlerCoupling<dim, spacedim> dof_handler_coupling(space_dh,
                                                                      dh);

  // The coupling mass matrix has shape (dh.n_dofs(), space_dh.n_dofs()),
  // which is the transpose of the layout used by the free functions in
  // non_matching/coupling.h.
  SparsityPattern sparsity;
  {
    DynamicSparsityPattern dsp(dh.n_dofs(), space_dh.n_dofs());
    dof_handler_coupling.create_coupling_mass_sparsity_pattern(quad, dsp);
    sparsity.copy_from(dsp);
  }
  SparseMatrix<double> coupling(sparsity);
  dof_handler_coupling.create_coupling_mass_matrix(quad, coupling);

  SparsityPattern mass_sparsity;
  {
    DynamicSparsityPattern dsp(dh.n_dofs(), dh.n_dofs());
    DoFTools::make_sparsity_pattern(dh, dsp);
    mass_sparsity.copy_from(dsp);
  }
  SparseMatrix<double> mass_matrix(mass_sparsity);
  MatrixTools::create_mass_matrix(dh, quad, mass_matrix);

  SparseDirectUMFPACK mass_matrix_inv;
  mass_matrix_inv.factorize(mass_matrix);

  // now take ones in space, project them onto the immersed space,
  // get back ones, and check for the error.
  Vector<double> space_ones(space_dh.n_dofs());
  Vector<double> ones(dh.n_dofs());

  space_ones = 1.0;
  coupling.vmult(ones, space_ones);
  mass_matrix_inv.solve(ones);

  Vector<double> real_ones(dh.n_dofs());
  real_ones = 1.0;
  ones -= real_ones;

  deallog << "Error on constants: " << ones.l2_norm() << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
