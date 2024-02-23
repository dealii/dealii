// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



/* Compare the kelly estimator for the same square but in codimension 0 and the
   other in codimension 1, Both should return the same estimator */



#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim, int spacedim>
void
check()
{
  Functions::CosineFunction<spacedim> function;
  Triangulation<dim, spacedim>        tr;

  GridGenerator::hyper_cube(tr, -1, 1);

  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();


  FE_Q<dim, spacedim>       element(QIterated<1>(QTrapezoid<1>(), 3));
  DoFHandler<dim, spacedim> dof(tr);
  dof.distribute_dofs(element);

  MappingQ<dim, spacedim> mapping(3);
  QGauss<dim - 1>         q_face(4);

  std::map<types::boundary_id, const Function<spacedim> *> neumann_bc;
  neumann_bc[0] = &function;

  Vector<double> v(dof.n_dofs());
  VectorTools::interpolate(mapping, dof, function, v);

  Vector<float> error(tr.n_active_cells());

  KellyErrorEstimator<dim, spacedim>::estimate(
    mapping, dof, q_face, neumann_bc, v, error);

  deallog << "Estimated error:" << std::endl;
  for (unsigned int i = 0; i < error.size(); ++i)
    deallog << error(i) * 100 << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  deallog.push("2d_2");
  check<2, 2>();
  deallog.pop();
  deallog.push("2d_3");
  check<2, 3>();
  deallog.pop();
}
