//----------------------------  denis_1.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002 by the deal.II authors and Anna Schneebeli
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  denis_1.cc  ---------------------------


// check for a bug in DerivativeApproximation::approximate_gradient
//
// this program is a modified version of one by
// Denis Danilov <danilovdenis@yandex.ru>, 


#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <numerics/data_out.h>
#include <fe/fe_q.h>
#include <numerics/derivative_approximation.h>
#include <fstream>
#include <iostream>

int main() 
{
  std::ofstream logfile("denis_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Triangulation<2>   triangulation;
  FE_Q<2>            fe(2);
  DoFHandler<2>      dof_handler(triangulation);
  Vector<double>       phi_solution;
  Vector<float>        gradient_phi;
  float                gradient_phi_min, gradient_phi_max;
  double delta = 0.05;

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(5);

  dof_handler.distribute_dofs(fe);

  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  
  phi_solution.reinit(dof_handler.n_dofs());

  for(; cell!=endc; ++cell)
    {
      for(unsigned int vertex=0;
  	   vertex < GeometryInfo<2>::vertices_per_cell;
  	   ++vertex)                               
  	{
  	  double x, y, r;
	  x = cell->vertex(vertex)(0);
  	  y = cell->vertex(vertex)(1);
	  r = sqrt( x*x + y*y );
	  phi_solution( cell->vertex_dof_index(vertex,0) )
	    = 0.5 * ( 1 - tanh( ( r - 0.5 )
				/( 2*M_SQRT2 * delta ) ) ) ;
  	}
    }

  gradient_phi.reinit(triangulation.n_active_cells());
  DerivativeApproximation::approximate_gradient(dof_handler, phi_solution, gradient_phi);

  gradient_phi_min = 1e30;
  gradient_phi_max = -1;

  cell = dof_handler.begin_active();  
  for(unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
    {
      if( gradient_phi(cell_no) < gradient_phi_min )
	gradient_phi_min =  gradient_phi(cell_no);
      
      if( gradient_phi(cell_no) > gradient_phi_max )
	gradient_phi_max =  gradient_phi(cell_no);
      
    }

  deallog << "gradient_phi_min: " << gradient_phi_min << std::endl;
  deallog << "gradient_phi_max: " << gradient_phi_max << std::endl;

  //////////////////////
  
  
};
