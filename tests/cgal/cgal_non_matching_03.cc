// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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

// Compute the coupling mass matrix <v_i,q_j> and check its correct by computing
// the measure of the embedded grid.

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/non_matching/coupling.h>

#include "../tests.h"
using namespace dealii;


template <int dim, int spacedim>
void
test()
{
  constexpr int                     degree = 3;
  constexpr double                  left   = 0.5;
  constexpr double                  right  = .75;
  Triangulation<spacedim, spacedim> space_tria;
  Triangulation<dim, spacedim>      embedded_tria;

  GridGenerator::hyper_cube(space_tria, 0., 1.);
  GridGenerator::hyper_cube(embedded_tria, left, right);
  space_tria.refine_global(2);
  embedded_tria.refine_global(2);

  DoFHandler<spacedim>      space_dh(space_tria);
  DoFHandler<dim, spacedim> embedded_dh(embedded_tria);

  FE_Q<spacedim>      fe_space(1);
  FE_Q<dim, spacedim> fe_embedded(1);

  space_dh.distribute_dofs(fe_space);
  embedded_dh.distribute_dofs(fe_embedded);


  auto space_cache =
    std::make_unique<GridTools::Cache<spacedim>>(space_tria); // Q1 mapping
  auto embedded_cache = std::make_unique<GridTools::Cache<dim, spacedim>>(
    embedded_tria); // Q1 mapping

  // Compute Quadrature formulas on the intersections of the two
  const auto vec_info =
    NonMatching::compute_intersection(*space_cache, *embedded_cache, degree);

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> coupling_matrix(sparsity_pattern);

  AffineConstraints<double> constraints;
  AffineConstraints<double> embedded_constraints;
  DynamicSparsityPattern    dsp(space_dh.n_dofs(), embedded_dh.n_dofs());
  NonMatching::create_coupling_sparsity_pattern_with_exact_intersections(
    vec_info,
    space_dh,
    embedded_dh,
    dsp,
    constraints,
    ComponentMask(),
    ComponentMask(),
    embedded_constraints);

  sparsity_pattern.copy_from(dsp);
  coupling_matrix.reinit(sparsity_pattern);

  NonMatching::create_coupling_mass_matrix_with_exact_intersections(
    space_dh,
    embedded_dh,
    vec_info,
    coupling_matrix,
    constraints,
    ComponentMask(),
    ComponentMask(),
    MappingQ1<3>(),
    MappingQ1<3>(),
    embedded_constraints);

  Vector<double> ones_space(space_dh.n_dofs());
  Vector<double> ones_embedded(embedded_dh.n_dofs());
  ones_space    = 1.0;
  ones_embedded = 1.0;
  const double result =
    coupling_matrix.matrix_scalar_product(ones_space, ones_embedded);
  deallog << "Result with coupling matrix: " << std::setprecision(10) << result
          << std::endl;

  deallog << "Expected : " << std::setprecision(10)
          << std::pow(right - left, spacedim) << std::endl;
}

int
main()
{
  initlog();



  test<3, 3>();
}
