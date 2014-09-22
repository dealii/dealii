// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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



/* Author: Markus Buerg, 2012 */
/* Purpose: Check FEFieldFunction for hp::DoFHandler. */



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/numerics/fe_field_function.h>

#include <fstream>



template <int dim>
void
check ()
{
  Triangulation<dim> triangulation;

  GridGenerator::subdivided_hyper_cube (triangulation, 2);
  hp::FECollection<dim> fe_collection;

  for (unsigned int i = 1; i <= triangulation.n_active_cells (); ++i)
    fe_collection.push_back (FE_Q<dim> (i));

  hp::DoFHandler<dim> dof_handler (triangulation);

  dof_handler.distribute_dofs (fe_collection);

  Vector<double> vector (dof_handler.n_dofs ());

  for (unsigned int i = 0; i < dof_handler.n_dofs (); ++i)
    vector (i) = i;

  Functions::FEFieldFunction<dim, hp::DoFHandler<dim> >
  fe_field (dof_handler, vector);
  QGauss<dim> quadrature (5);

  deallog << "values:" <<std::endl;

  std::vector<double> values (quadrature.size ());

  fe_field.value_list (quadrature.get_points (), values);

  for (unsigned int q_point = 0; q_point < quadrature.size (); ++q_point)
    deallog << values[q_point] << std::endl;

  deallog << "gradients:" <<std::endl;

  std::vector<Tensor<1, dim> > gradients (quadrature.size ());

  fe_field.gradient_list (quadrature.get_points (), gradients);

  for (unsigned int q_point = 0; q_point < quadrature.size (); ++q_point)
    deallog << gradients[q_point] << std::endl;
}


int main ()
{
  initlog();
  deallog << std::setprecision (2);
  deallog << std::fixed;

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
