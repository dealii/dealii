// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that we pick the correct mapping in KellyErrorEstimator
// for hp applications.


#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // We investigate the effect of different mappings in an hp-context.
  // Thus we only feed the MappingCollection with different mappings,
  // and keep the FEs and quadratures the same.
  hp::FECollection<dim>      fes(FE_Q<dim>(1), FE_Q<dim>(1));
  hp::QCollection<dim - 1>   quads_face(QGauss<dim - 1>(2), QGauss<dim - 1>(2));
  hp::MappingCollection<dim> mappings(MappingQ<dim>(1), MappingQ<dim>(2));

  // Set up hyper_ball and assign different active FE indices
  // to inner and outer cells.
  Triangulation<dim> tria;
  DoFHandler<dim>    dofh(tria);

  GridGenerator::hyper_ball(tria);
  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->center() != Point<dim>())
      cell->set_active_fe_index(1);
  tria.refine_global(1);

  dofh.distribute_dofs(fes);

  // Interpolate a transcendental function and estimate the error via Kelly
  // with different mappings.
  Functions::CosineFunction<dim> function;
  Vector<double>                 interpolation(dofh.n_dofs());

  Vector<float> error(tria.n_active_cells());

  // MappingQ1
  hp::MappingCollection<dim> mapping_q1(mappings[0]);
  VectorTools::interpolate(mapping_q1, dofh, function, interpolation);
  KellyErrorEstimator<dim>::estimate(
    mapping_q1, dofh, quads_face, {}, interpolation, error);
  deallog << "MappingQ1  : global error estimate=" << error.l2_norm()
          << std::endl;

  // MappingQ2
  hp::MappingCollection<dim> mapping_q2(mappings[1]);
  VectorTools::interpolate(mapping_q2, dofh, function, interpolation);
  KellyErrorEstimator<dim>::estimate(
    mapping_q2, dofh, quads_face, {}, interpolation, error);
  deallog << "MappingQ2  : global error estimate=" << error.l2_norm()
          << std::endl;

  // MappingQ1&2
  VectorTools::interpolate(mappings, dofh, function, interpolation);
  KellyErrorEstimator<dim>::estimate(
    mappings, dofh, quads_face, {}, interpolation, error);
  deallog << "MappingQ1&2: global error estimate=" << error.l2_norm()
          << std::endl;
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
