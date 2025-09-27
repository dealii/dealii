// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Ensure that FEFaceEvaluation can be collected in a vector

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include "../tests.h"


int
main()
{
  initlog();

  constexpr int      dim = 2;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingQ<dim> mapping(1);

  MatrixFree<dim, double>                 matrix_free;
  MatrixFree<dim, double>::AdditionalData mf_data;
  mf_data.mapping_update_flags_boundary_faces = update_values;
  matrix_free.reinit(
    mapping, dof_handler, AffineConstraints<double>(), QGauss<1>(2), mf_data);
  FEFaceEvaluation<dim, -1, 0, 1, double>              evaluator(matrix_free);
  std::vector<FEFaceEvaluation<dim, -1, 0, 1, double>> evaluators(2, evaluator);
  evaluators[0].reinit(0);
  evaluators[1].reinit(0);

  deallog << "OK" << std::endl;
}
