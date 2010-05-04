//----------------------------  crash_19.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  crash_19.cc  ---------------------------


// VectorTools::interpolate_boundary_values produced an exception when
// used with hp::DoFHandler in 1d. Test that this is no longer the case


#include "../tests.h"
#include <base/logstream.h>
#include <fstream>
std::ofstream logfile("crash_19/output");


#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/table_handler.h>
#include <base/thread_management.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_refinement.h>
#include <dofs/dof_handler.h>
#include <lac/constraint_matrix.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <hp/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

#include <iostream>
#include <fstream>
#include <list>
#include <sstream>


template <int dim>
class ExactSolution: public Function<dim> {
  public:
    ExactSolution () {}
    virtual double value (const Point<dim>& p, const unsigned int) const
      {
	return p (0);
      }
 };


template <int dim>
void test ()
{
  Triangulation<dim>     triangulation;
  hp::FECollection<dim>      fe;
  fe.push_back (FE_Q<dim> (1));

  hp::DoFHandler<dim>        dof_handler (triangulation);

  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (2);
  deallog << "Number of active cells: "
	  << triangulation.n_active_cells()
	  << std::endl;
  deallog << "Total number of cells: "
	  << triangulation.n_cells()
	  << std::endl;

  dof_handler.distribute_dofs (fe);
  deallog << "Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << std::endl;

  ExactSolution<dim> exact_solution;
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            exact_solution,
                                            boundary_values);
  if (dim == 1)
    VectorTools::interpolate_boundary_values (dof_handler,
					      1,
					      exact_solution,
					      boundary_values);

  for (std::map<unsigned int,double>::iterator i=boundary_values.begin();
       i != boundary_values.end(); ++i)
    deallog << i->first << ' ' << i->second << std::endl;
}


int main ()
{
  try
    {
      logfile.precision(2);
      deallog << std::setprecision(2);

      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test<1> ();
      test<2> ();
      test<3> ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };

  return 0;
}
