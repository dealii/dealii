//---------------------------  create_point_source.cc  --------------------------
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
//---------------------------  create_point_source.cc  --------------------------


// Test the function VectorTools::create_point_source_vector.



#include "../tests.h"
#include<deal.II/dofs/dof_handler.h>
#include<deal.II/fe/fe_q.h>
#include<deal.II/fe/fe_system.h>
#include<deal.II/grid/grid_generator.h>
#include<deal.II/grid/tria.h>
#include<deal.II/numerics/vector_tools.h>

#include<fstream>




template <int dim>
void
check ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  FESystem<dim> fe (FE_Q<dim> (2), dim);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs (fe);
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
  std::ofstream logfile ("create_point_source/output");
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
