// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/non_matching/coupling.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


// Test that a coupling matrix can be constructed for each pair of dimension
// and immersed dimension, and check that quadratic functions are correctly
// projected.

template <int dim, int spacedim>
void
test()
{
  deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl;

  Triangulation<dim, spacedim>      tria;
  Triangulation<spacedim, spacedim> space_tria;

  GridGenerator::hyper_cube(tria, -.4, .3);
  GridGenerator::hyper_cube(space_tria, -1, 1);

  tria.refine_global(1);
  space_tria.refine_global(2);

  FE_Q<dim, spacedim>      fe(2);
  FE_Q<spacedim, spacedim> space_fe(2);

  deallog << "FE      : " << fe.get_name() << std::endl
          << "Space FE: " << space_fe.get_name() << std::endl;

  DoFHandler<dim, spacedim>      dh(tria);
  DoFHandler<spacedim, spacedim> space_dh(space_tria);

  dh.distribute_dofs(fe);
  space_dh.distribute_dofs(space_fe);

  deallog << "Dofs      : " << dh.n_dofs() << std::endl
          << "Space dofs: " << space_dh.n_dofs() << std::endl;

  QGauss<dim> quad(3); // Quadrature for coupling


  SparsityPattern sparsity;
  {
    DynamicSparsityPattern dsp(space_dh.n_dofs(), dh.n_dofs());
    NonMatching::create_coupling_sparsity_pattern(space_dh, dh, quad, dsp);
    sparsity.copy_from(dsp);
  }
  SparseMatrix<double> coupling(sparsity);
  NonMatching::create_coupling_mass_matrix(
    space_dh, dh, quad, coupling, AffineConstraints<double>());

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

  // now take the square function in space, project them onto the immersed
  // space, get back ones, and check for the error.
  Vector<double> space_square(space_dh.n_dofs());
  Vector<double> squares(dh.n_dofs());
  Vector<double> projected_squares(dh.n_dofs());

  VectorTools::interpolate(space_dh,
                           Functions::SquareFunction<spacedim>(),
                           space_square);
  VectorTools::interpolate(dh, Functions::SquareFunction<spacedim>(), squares);

  coupling.Tvmult(projected_squares, space_square);
  mass_matrix_inv.solve(projected_squares);

  projected_squares -= squares;

  deallog << "Error on squares: " << projected_squares.l2_norm() << std::endl;
}



int
main()
{
  initlog();
  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
