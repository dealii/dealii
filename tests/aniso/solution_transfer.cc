// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>
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
  std::ofstream logfile("output");
  deallog << std::setprecision (4);
  logfile << std::setprecision (4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  transfer<2>(logfile);
  transfer<3>(logfile);
}



