//----------------------------  crash_12.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  crash_12.cc  ---------------------------


// check the complex case described in the hp paper by playing through all
// sorts of arrangements of finite elements on one coarse and one refined cell
//
// this code in particular tests some compensating code in
// dof_tools.cc, where we have to make sure that we select a suitable
// set of master dofs. this is mostly trivial in 2d and for most fe
// combinations in 3d as well. the exceptions are that it doesn't work
// as easily in 3d for the combinations Q4/Q3, Q5/Q3, and
// Q5/Q4. Higher order finite elements in 3d will probably only
// exacerbate the problem, but the code there appears to be robust.

char logname[] = "crash_12/output";


#include "../tests.h"
#include <base/function.h>
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>

#include <grid/tria.h>
#include <hp/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <numerics/vectors.h>
#include <fe/fe_q.h>

#include <fstream>
#include <vector>



template <int dim>
void test ()
{
				   // create a mesh like this (viewed
				   // from top, if in 3d):
				   // *---*---*
				   // | 0 | 1 |
				   // *---*---*
  Triangulation<dim>     triangulation;
  std::vector<unsigned int> subdivisions (dim, 1);
  subdivisions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions,
                                             Point<dim>(),
					     (dim == 3 ?
					      Point<dim>(2,1,1) :
					      Point<dim>(2,1)));
  (++triangulation.begin_active())->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();
  
  hp::FECollection<dim> fe;
  fe.push_back (FE_Q<dim>(1));
  fe.push_back (FE_Q<dim>(2));
  fe.push_back (FE_Q<dim>(3));
  fe.push_back (FE_Q<dim>(4));
  fe.push_back (FE_Q<dim>(5));

  hp::DoFHandler<dim>        dof_handler(triangulation);

  for (unsigned int i=0; i<fe.size(); ++i)
    for (unsigned int j=0; j<fe.size(); ++j)
      {
	deallog << "Testing " << fe[i].get_name()
		<< " vs. " << fe[j].get_name()
		<< std::endl;
	
					 // set fe on coarse cell to 'i', on
					 // all fine cells to 'j'
	typename hp::DoFHandler<dim>::active_cell_iterator
	  cell = dof_handler.begin_active();
	cell->set_active_fe_index (i);
	++cell;

	for (; cell != dof_handler.end(); ++cell)
	  cell->set_active_fe_index (j);

	dof_handler.distribute_dofs (fe);
  
	ConstraintMatrix constraints;
	DoFTools::make_hanging_node_constraints (dof_handler,
						 constraints);
	constraints.close ();
	
	constraints.print (deallog.get_file_stream());
      }
}



int main ()
{
  std::ofstream logfile(logname);
  logfile.precision (4);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);


  test<2> ();
  test<3> ();
}

