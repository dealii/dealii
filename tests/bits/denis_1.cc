//----------------------------  denis_1.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  denis_1.cc  ---------------------------


// check for something in
// DerivativeApproximation::approximate_gradient. the original report
// stated it was a bug, but it was not (see the archives), but since
// the program is already there let's make use of it.
//
// this program is a modified version of one by
// Denis Danilov <danilovdenis@yandex.ru>, 


#include "../tests.h"
#include <base/logstream.h>
#include <base/function.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <numerics/vectors.h>
#include <fe/fe_q.h>
#include <numerics/derivative_approximation.h>
#include <fstream>
#include <iostream>


class F : public Function<2>
{
  public:
    virtual double value (const Point<2> &p, const unsigned int) const
      {
        double delta = 0.05;
        double x, y, r;
        x = p(0);
        y = p(1);
        r = sqrt( x*x + y*y );
        return 0.5 * ( 1 - tanh( ( r - 0.5 )
                                 /( 2*M_SQRT2 * delta ) ) ) ;
      };
};


int main() 
{
  std::ofstream logfile("denis_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  logfile.precision (2);
  logfile.setf(std::ios::fixed);  

  Triangulation<2>   triangulation;
  FE_Q<2>            fe(2);
  DoFHandler<2>      dof_handler(triangulation);
  Vector<double>       phi_solution;
  Vector<float>        gradient_phi;
  float                gradient_phi_min, gradient_phi_max;

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(5);

  dof_handler.distribute_dofs(fe);

  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  
  phi_solution.reinit(dof_handler.n_dofs());
  VectorTools::interpolate (dof_handler, F(), phi_solution);

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
}
