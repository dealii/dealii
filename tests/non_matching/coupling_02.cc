// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#include "../tests.h"


// Test that a coupling matrix can be constructed for each pair of dimension and
// immersed dimension, with a non trivial component selection in the space_dh,
// and check that constants are projected correctly.

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

  FE_Q<dim, spacedim>          fe(1);
  FESystem<spacedim, spacedim> space_fe(FE_Q<spacedim, spacedim>(1),
                                        spacedim + 1);
  ComponentMask                space_mask(spacedim + 1, false);
  space_mask.set(spacedim, true);

  deallog << "FE      : " << fe.get_name() << std::endl
          << "Space FE: " << space_fe.get_name() << std::endl;

  DoFHandler<dim, spacedim>      dh(tria);
  DoFHandler<spacedim, spacedim> space_dh(space_tria);

  dh.distribute_dofs(fe);
  space_dh.distribute_dofs(space_fe);

  deallog << "Dofs      : " << dh.n_dofs() << std::endl
          << "Space dofs: " << space_dh.n_dofs() << std::endl;

  QGauss<dim> quad(3); // Quadrature for coupling

  AffineConstraints<double> constraints;

  SparsityPattern sparsity;
  {
    DynamicSparsityPattern dsp(space_dh.n_dofs(), dh.n_dofs());
    NonMatching::create_coupling_sparsity_pattern(
      space_dh, dh, quad, dsp, constraints, space_mask);
    sparsity.copy_from(dsp);
  }
  SparseMatrix<double> coupling(sparsity);
  NonMatching::create_coupling_mass_matrix(
    space_dh, dh, quad, coupling, constraints, space_mask);

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
  coupling.Tvmult(ones, space_ones);
  mass_matrix_inv.solve(ones);

  Vector<double> real_ones(dh.n_dofs());
  real_ones = 1.0;
  ones -= real_ones;

  deallog << "Error on constants: " << ones.l2_norm() << std::endl;
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
