// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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



/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/derivative_approximation.h>

#include <fstream>



template <int dim>
void
check ()
{
  Functions::CosineFunction<dim> cosine;

  Triangulation<dim> tr;
  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  FE_Q<dim> element(QIterated<1>(QTrapez<1>(),3));
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
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);

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
