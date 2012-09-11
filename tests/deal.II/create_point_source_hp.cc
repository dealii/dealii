//--------------------------  create_point_source_hp.cc  -------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------  create_point_source_hp.cc  -------------------------


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
  std::ofstream logfile ("create_point_source_hp/output");
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
