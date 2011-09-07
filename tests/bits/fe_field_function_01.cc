//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// Test the FEFieldFunction class

#include "../tests.h"
#include <fstream>

// all include files you need here
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vectors.h>

double abs_zero(double a)
{
  if( std::abs(a) < 1e-10)
    return 0;
  else
    return a;
}

template <int dim>
void test() {
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(8/dim);

  FE_Q<dim> fe(1);
  DoFHandler<dim> dh(tria);

  dh.distribute_dofs(fe);
  deallog << "Ndofs :" << dh.n_dofs() << std::endl;

  Functions::CosineFunction<dim> ff;

  Vector<double> v1(dh.n_dofs()), v2(dh.n_dofs());

  VectorTools::interpolate(dh, ff, v1);
  deallog << "V norm: " << v1.l2_norm() << std::endl;

  Functions::FEFieldFunction<dim, DoFHandler<dim>, Vector<double> >
    fef(dh, v1);

  VectorTools::interpolate(dh, fef, v2);

  v2.add(-1, v1);
  deallog << "Interpolation error: " << abs_zero(v2.l2_norm())
	  << std::endl;

  Vector<double> error(tria.n_active_cells());
  QGauss<dim> quad(2);
  VectorTools::integrate_difference(dh, v1, ff,
				    error, quad,
				    VectorTools::H1_norm);
  deallog << "H1 Interpolation error: "
	  << abs_zero(error.l2_norm()) << std::endl;
  error = 0;
  VectorTools::integrate_difference(dh, v1, fef,
				    error, quad,
				    VectorTools::H1_norm);
  deallog << "H1 Interpolation error with fef: "
	  << abs_zero(error.l2_norm()) << std::endl;

}

int main ()
{
  std::ofstream logfile("fe_field_function_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}

