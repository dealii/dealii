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
#include <deal.II/fe/fe_q.h>
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
    double f=sin(p[0]*4);
    if (dim>1)
      f*=cos(p[1]*4);
    if (dim>2)
      f*=exp(p[2]*4);
    return f;
  };
};


template <int dim>
void transfer(std::ostream &out)
{
  MyFunction<dim> function;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5-dim);
  FE_Q<dim> fe_q(1);
  FE_DGQ<dim> fe_dgq(1);
  DoFHandler<dim> q_dof_handler(tria);
  DoFHandler<dim> dgq_dof_handler(tria);
  Vector<double> q_solution;
  Vector<double> dgq_solution;
  MappingQ1<dim> mapping;
  DataOut<dim> q_data_out, dgq_data_out;
  ConstraintMatrix cm;
  cm.close();

  q_dof_handler.distribute_dofs (fe_q);
  q_solution.reinit(q_dof_handler.n_dofs());

  dgq_dof_handler.distribute_dofs (fe_dgq);
  dgq_solution.reinit(dgq_dof_handler.n_dofs());

  VectorTools::interpolate (mapping, q_dof_handler, function, q_solution);
  VectorTools::project (mapping, dgq_dof_handler, cm,
                        QGauss<dim>(3), function, dgq_solution);

  q_data_out.attach_dof_handler (q_dof_handler);
  q_data_out.add_data_vector (q_solution, "solution");
  q_data_out.build_patches ();
  deallog << "Initial solution, FE_Q" << std::endl << std::endl;
  q_data_out.write_gnuplot (out);

  dgq_data_out.attach_dof_handler (dgq_dof_handler);
  dgq_data_out.add_data_vector (dgq_solution, "solution");
  dgq_data_out.build_patches ();
  deallog << "Initial solution, FE_DGQ" << std::endl << std::endl;
  dgq_data_out.write_gnuplot (out);

  SolutionTransfer<dim> q_soltrans(q_dof_handler);
  SolutionTransfer<dim> dgq_soltrans(dgq_dof_handler);

  // test a): pure refinement
  typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active(),
                                                    endc=tria.end();
  ++cell;
  ++cell;
  for (; cell!=endc; ++cell)
    cell->set_refine_flag();

  tria.prepare_coarsening_and_refinement();
  q_soltrans.prepare_for_pure_refinement();
  dgq_soltrans.prepare_for_pure_refinement();
  tria.execute_coarsening_and_refinement();
  q_dof_handler.distribute_dofs (fe_q);
  dgq_dof_handler.distribute_dofs (fe_dgq);

  Vector<double> tmp_q (q_dof_handler.n_dofs());
  q_soltrans.refine_interpolate(q_solution, tmp_q);
  q_solution.reinit (q_dof_handler.n_dofs());
  q_solution = tmp_q;

  Vector<double> tmp_dgq (dgq_dof_handler.n_dofs());
  dgq_soltrans.refine_interpolate(dgq_solution, tmp_dgq);
  dgq_solution.reinit (dgq_dof_handler.n_dofs());
  dgq_solution = tmp_dgq;


  q_data_out.clear_data_vectors();
  q_data_out.add_data_vector (q_solution, "solution");
  q_data_out.build_patches ();
  deallog << "Interpolated/tranferred solution after pure refinement, FE_Q"
          << std::endl << std::endl;
  q_data_out.write_gnuplot (out);

  dgq_data_out.clear_data_vectors();
  dgq_data_out.add_data_vector (dgq_solution, "solution");
  dgq_data_out.build_patches ();
  deallog << "Interpolated/tranferred solution after pure refinement, FE_DGQ"
          << std::endl << std::endl;
  dgq_data_out.write_gnuplot (out);

  // test b): with coarsening
  q_soltrans.clear();
  dgq_soltrans.clear();

  cell=tria.begin_active(tria.n_levels()-1);
  endc=tria.end(tria.n_levels()-1);
  cell->set_refine_flag();
  ++cell;
  for (; cell!=endc; ++cell)
    cell->set_coarsen_flag();
  Vector<double> q_old_solution=q_solution,
                 dgq_old_solution=dgq_solution;
  tria.prepare_coarsening_and_refinement();
  q_soltrans.prepare_for_coarsening_and_refinement(q_old_solution);
  dgq_soltrans.prepare_for_coarsening_and_refinement(dgq_old_solution);
  tria.execute_coarsening_and_refinement();
  q_dof_handler.distribute_dofs (fe_q);
  dgq_dof_handler.distribute_dofs (fe_dgq);
  q_solution.reinit(q_dof_handler.n_dofs());
  dgq_solution.reinit(dgq_dof_handler.n_dofs());
  q_soltrans.interpolate(q_old_solution, q_solution);
  dgq_soltrans.interpolate(dgq_old_solution, dgq_solution);

  q_data_out.clear_data_vectors();
  q_data_out.add_data_vector (q_solution, "solution");
  q_data_out.build_patches ();
  deallog << "Interpolated/tranferred solution after coarsening and refinement, FE_Q"
          << std::endl << std::endl;
  q_data_out.write_gnuplot (out);

  dgq_data_out.clear_data_vectors();
  dgq_data_out.add_data_vector (dgq_solution, "solution");
  dgq_data_out.build_patches ();
  deallog << "Interpolated/tranferred solution after coarsening and refinement, FE_DGQ"
          << std::endl << std::endl;
  dgq_data_out.write_gnuplot (out);
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << "   1D solution transfer" << std::endl;
  transfer<1>(logfile);

  deallog << "   2D solution transfer" << std::endl;
  transfer<2>(logfile);

  deallog << "   3D solution transfer" << std::endl;
  transfer<3>(logfile);
}



