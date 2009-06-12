//----------------------------  solution_transfer.cc  ---------------------------
//    fe_data_test.cc,v 1.14 2003/11/28 11:52:35 guido Exp
//    Version: 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solution_transfer.cc  ---------------------------


#include "../tests.h"
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>

#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <numerics/vectors.h>
#include <numerics/solution_transfer.h>
#include <numerics/data_out.h>
#include <fe/fe_dgq.h>
#include <fe/mapping_q1.h>
#include <fstream>
#include <iostream>
#include <vector>


template<int dim>
class MyFunction : public Function<dim>
{
  public:
    MyFunction () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int) const
      {
	double ret_value=sin(p[0]*4)*cos(p[1]*4);
	if (dim==3)
	  ret_value*=sin(5*p[2]+1);
	return ret_value;
      };
};


template <int dim>
void transfer(std::ostream &out)
{
  MyFunction<dim> function;
  Triangulation<dim> tria(Triangulation<dim>::allow_anisotropic_smoothing);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  FE_DGQ<dim> fe(1);
  DoFHandler<dim> dof_handler(tria);
  Vector<double> solution;
  MappingQ1<dim> mapping;
  DataOut<dim> data_out;
  
  dof_handler.distribute_dofs (fe);
  solution.reinit(dof_handler.n_dofs());

  VectorTools::interpolate (mapping, dof_handler, function, solution);

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  deallog << "Initial solution" << std::endl << std::endl;
  data_out.write_gnuplot (out);

  SolutionTransfer<dim> soltrans(dof_handler);

				   // test a): pure refinement
  typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active(),
						    endc=tria.end();
  for (; cell!=endc; ++cell)
    cell->set_refine_flag(RefinementCase<dim>::cut_x);

  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_pure_refinement();
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs (fe);

  Vector<double> new_solution (dof_handler.n_dofs());
  soltrans.refine_interpolate(solution, new_solution);
  solution.reinit (dof_handler.n_dofs());
  solution = new_solution;

  data_out.clear_data_vectors();
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  deallog << "Interpolated/tranferred solution after pure refinement" << std::endl << std::endl;
  data_out.write_gnuplot (out);

				   // test b): with coarsening
  SolutionTransfer<dim> soltrans2(dof_handler);
  cell=tria.begin_active(tria.n_levels()-1);
  endc=tria.end(tria.n_levels()-1);
  for (; cell!=endc; ++cell)
    cell->set_coarsen_flag();
  Vector<double> old_solution=solution;
  tria.prepare_coarsening_and_refinement();
  soltrans2.prepare_for_coarsening_and_refinement(old_solution);
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs (fe);
  solution.reinit(dof_handler.n_dofs());
  soltrans2.interpolate(old_solution, solution);

  data_out.clear_data_vectors();
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  deallog << "Interpolated/tranferred solution after coarsening" << std::endl << std::endl;
  data_out.write_gnuplot (out);

}


int main()
{
  std::ofstream logfile("solution_transfer/output");
  deallog << std::setprecision (4);
  logfile << std::setprecision (4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  transfer<2>(logfile);
  transfer<3>(logfile);
}



