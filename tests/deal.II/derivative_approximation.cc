//----------------------------  derivative_approximation.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  derivative_approximation.cc  ---------------------------


/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */



#include <base/logstream.h>
#include <base/function_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <numerics/vectors.h>
#include <numerics/derivative_approximation.h>

#include <fstream>



template <int dim>
void
check ()
{
  CosineFunction<dim> cosine;
  
  Triangulation<dim> tr;  
  if (dim==2)
    GridGenerator::hyper_ball(tr);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);
  
  FE_Q<dim> element(3);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  Vector<double> v (dof.n_dofs());
  VectorTools::interpolate (dof, cosine, v);

  Vector<float> gradient (tr.n_active_cells());
  Vector<float> second (tr.n_active_cells());

  DerivativeApproximation::
    approximate_gradient (dof, v, gradient);
  DerivativeApproximation::
    approximate_second_derivative (dof, v, second);
  
  deallog << "Approximated gradient:" << std::endl;
  for (unsigned int i=0; i<gradient.size(); ++i)
    deallog << gradient(i)*100 << std::endl;
  
  deallog << "Approximated second derivative:" << std::endl;
  for (unsigned int i=0; i<gradient.size(); ++i)
    deallog << second(i)*100 << std::endl;
}

int main ()
{
  std::ofstream logfile ("derivative_approximation.output");
  logfile.precision (2);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
