//----------------------------  vectors_point_source_hp.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vectors_point_source_hp.cc  ---------------------------


// like deal.II/vectors_point_source_hp, but for hp objects. here, each hp object has only a
// single component, so we expect exactly the same output as for the old test.
// vectors_point_source_hp_hp tests for different finite elements


#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function_lib.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <hp/dof_handler.h>
#include <dofs/dof_tools.h>
#include <lac/constraint_matrix.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <hp/fe_collection.h>
#include <hp/q_collection.h>
#include <fe/mapping_q.h>
#include <hp/mapping_collection.h>
#include <numerics/vectors.h>

#include <fstream>


template <int dim>
void
check ()
{
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

				   // create a system element composed
				   // of one Q1 and one Q2 element
  hp::FECollection<dim> element;
  for (unsigned int i=1; i<7-dim; ++i)
    element.push_back (FE_Q<dim>(i));
  hp::DoFHandler<dim> dof(tr);
  for (typename hp::DoFHandler<dim>::active_cell_iterator
	 cell = dof.begin_active(); cell!=dof.end(); ++cell)
    cell->set_active_fe_index (rand() % element.size());
  
  dof.distribute_dofs(element);

				   // use a more complicated mapping
				   // of the domain and a quadrature
				   // formula suited to the elements
				   // we have here
  hp::MappingCollection<dim> mapping;
  for (unsigned int i=1; i<7-dim; ++i)
    mapping.push_back (MappingQ<dim>(i+1));

  Vector<double> rhs (dof.n_dofs());
  VectorTools::create_point_source_vector (mapping, dof,
					   tr.begin()->center(),
					   rhs);
  for (unsigned int i=0; i<rhs.size(); ++i)
    deallog << rhs(i) << std::endl;
}



int main ()
{
  std::ofstream logfile ("vectors_point_source_hp_01/output");
  logfile.precision (4);
  logfile.setf(std::ios::fixed);  
  deallog << std::setprecision(4) << std::fixed;
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
