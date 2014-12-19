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



// Test the function VectorTools::create_point_source_vector for hp.



#include "../tests.h"
#include<deal.II/fe/fe_q.h>
#include<deal.II/fe/fe_system.h>
#include<deal.II/grid/grid_generator.h>
#include<deal.II/grid/tria.h>
#include<deal.II/hp/dof_handler.h>
#include<deal.II/hp/fe_collection.h>
#include<deal.II/numerics/vector_tools.h>

#include<fstream>




template <int dim>
void
check ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  hp::FECollection<dim> fe_collection;

  for (unsigned int i = 1; i <= tria.n_active_cells (); ++i)
    fe_collection.push_back (FESystem<dim> (FE_Q<dim> (i), dim));

  hp::DoFHandler<dim> dof (tria);
  dof.distribute_dofs (fe_collection);
  Point<dim> orientation;
  Point<dim> p (tria.begin_active ()->center ());

  for (unsigned int i = 0; i < dim; ++i)
    orientation (i) = i;

  Vector<double> vector (dof.n_dofs ());

  VectorTools::create_point_source_vector (dof, p, orientation, vector);

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    deallog << vector (i) << std::endl;
}



int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
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
