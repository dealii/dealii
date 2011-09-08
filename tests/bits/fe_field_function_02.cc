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


// Test the FEFieldFunction class. This class was not thread-safe at one point
// because it keeps a cache on the side that was invalidated when different
// threads kept pouncing on it. Patrick Sodre wrote a fix for that that's
// based on a thread-local variable instead of the single cache.

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

				   // project the discrete function fef back
				   // onto the finite element space. this
				   // should be the identity operation and
				   // consequently subtracting the vector from
				   // itself should yield zero
  {
    ConstraintMatrix cm;
    cm.close();
    VectorTools::project(dh, cm, QGauss<dim>(3), fef, v2);
  }
  v2.add(-1, v1);

  deallog << "Projection error: " << abs_zero(v2.l2_norm())
	  << std::endl;
  Assert (v2.l2_norm() < 1e-12, ExcInternalError());
}

int main ()
{
  std::ofstream logfile("fe_field_function_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}

