//----------------------------  mesh_3d_12.cc  ---------------------------
//    mesh_3d_12.cc,v 1.2 2003/10/16 14:16:40 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_12.cc  ---------------------------


// the Kelly estimator used to fall over when presented with 3d meshes
// with mis-oriented faces. the actual assertion where we fail is
// again tested isolated in mesh_3d_13 (we limit the number of
// refinements here to reduce run time; it is not limited there to the
// same degree)

#include "../tests.h"
#include "mesh_3d.h"

#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <base/function_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <numerics/error_estimator.h>
#include <numerics/vectors.h>

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
                                    QGauss2<2>(),
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
  std::ofstream logfile("mesh_3d_12.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

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

  
  
