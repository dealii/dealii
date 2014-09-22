// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// the Kelly estimator used to fall over when presented with 3d meshes
// with mis-oriented faces. the actual assertion where we fail is
// again tested isolated in mesh_3d_13 (we limit the number of
// refinements here to reduce run time; it is not limited there to the
// same degree)

#include "../tests.h"
#include "mesh_3d.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>




void check_this (Triangulation<3> &tria)
{
  FE_Q<3> fe(1);
  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  Vector<double> u(dof_handler.n_dofs());
  Vector<float>  e(tria.n_active_cells());

  VectorTools::interpolate (dof_handler,
                            Functions::SquareFunction<3>(),
                            u);

  KellyErrorEstimator<3>::estimate (dof_handler,
                                    QGauss<2>(2),
                                    FunctionMap<3>::type(),
                                    u,
                                    e);

  deallog << "  " << e.l1_norm() << std::endl;
  deallog << "  " << e.l2_norm() << std::endl;
  deallog << "  " << e.linfty_norm() << std::endl;
}


void check (Triangulation<3> &tria)
{
  (++tria.begin_active())->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();

  deallog << "Initial check" << std::endl;
  check_this (tria);

  for (unsigned int r=0; r<2; ++r)
    {
      tria.refine_global (1);
      deallog << "Check " << r << std::endl;
      check_this (tria);
    }
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {
    Triangulation<3> coarse_grid;
    create_two_cubes (coarse_grid);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_L_shape (coarse_grid);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    check (coarse_grid);
  }

}



