// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// continuous projection of a function on the surface of a hypersphere

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <string>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  FE_Q<dim, spacedim>       fe(1);
  DoFHandler<dim, spacedim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);
  Functions::CosineFunction<spacedim> cosine;
  Vector<double>                      x(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, cosine, x);

  const double mean =
    VectorTools::compute_mean_value(dof_handler, QGauss<dim>(2), x, 0);
  // we have a symmetric domain and a symmetric function. the result
  // should be close to zero
  AssertThrow(std::fabs(mean) < 1e-15, ExcInternalError());
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test<1, 2>();
  test<2, 3>();
}
