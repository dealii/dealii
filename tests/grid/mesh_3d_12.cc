// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// the Kelly estimator used to fall over when presented with 3d meshes
// with mis-oriented faces. the actual assertion where we fail is
// again tested isolated in mesh_3d_13 (we limit the number of
// refinements here to reduce run time; it is not limited there to the
// same degree)

#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "mesh_3d.h"



void
check_this(Triangulation<3> &tria)
{
  FE_Q<3>       fe(1);
  DoFHandler<3> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> u(dof_handler.n_dofs());
  Vector<float>  e(tria.n_active_cells());

  VectorTools::interpolate(dof_handler, Functions::SquareFunction<3>(), u);

  KellyErrorEstimator<3>::estimate(
    dof_handler,
    QGauss<2>(2),
    std::map<types::boundary_id, const Function<3> *>(),
    u,
    e);

  deallog << "  " << static_cast<double>(e.l1_norm()) << std::endl;
  deallog << "  " << static_cast<double>(e.l2_norm()) << std::endl;
  deallog << "  " << static_cast<double>(e.linfty_norm()) << std::endl;
}


void
check(Triangulation<3> &tria)
{
  (std::next(tria.begin_active()))->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  deallog << "Initial check" << std::endl;
  check_this(tria);

  for (unsigned int r = 0; r < 2; ++r)
    {
      tria.refine_global(1);
      deallog << "Check " << r << std::endl;
      check_this(tria);
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(8);

  {
    Triangulation<3> coarse_grid;
    create_two_cubes(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_L_shape(coarse_grid);
    check(coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball(coarse_grid);
    coarse_grid.reset_manifold(0);
    check(coarse_grid);
  }
}
