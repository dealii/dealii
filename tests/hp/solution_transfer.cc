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
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <fstream>
#include <iostream>
#include <vector>

// a linear function that should be transferred exactly with Q1 and Q2
// elements
template<int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction () : Function<dim>() {};

  virtual double value (const Point<dim>   &p,
                        const unsigned int) const
  {
    double f=0.25 + 2 * p[0];
    if (dim>1)
      f+=0.772 * p[1];
    if (dim>2)
      f-=3.112 * p[2];
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
  const unsigned int max_degree = 6-dim;
  hp::FECollection<dim> fe_q;
  hp::FECollection<dim> fe_dgq;
  for (unsigned int deg=1; deg<=max_degree; ++deg)
    {
      fe_q.push_back (FE_Q<dim>(deg));
      fe_dgq.push_back (FE_DGQ<dim>(deg));
    }
  hp::DoFHandler<dim> q_dof_handler(tria);
  hp::DoFHandler<dim> dgq_dof_handler(tria);
  Vector<double> q_solution;
  Vector<double> dgq_solution;
  MappingQ1<dim> mapping;

  // refine a few cells
  typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active(),
                                                    endc=tria.end();
  ++cell;
  ++cell;
  for (; cell!=endc; ++cell)
    cell->set_refine_flag();
  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();

  // randomly assign FE orders
  unsigned int counter = 0;
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell=q_dof_handler.begin_active(),
    celldg=dgq_dof_handler.begin_active(),
    endc=q_dof_handler.end();
    for (; cell!=endc; ++cell, ++celldg, ++counter)
      {
        if (counter < 15)
          cell->set_active_fe_index(1);
        else
          cell->set_active_fe_index(Testing::rand()%max_degree);
        if (counter < 15)
          celldg->set_active_fe_index(1);
        else
          celldg->set_active_fe_index(Testing::rand()%max_degree);
      }
  }

  q_dof_handler.distribute_dofs (fe_q);
  q_solution.reinit(q_dof_handler.n_dofs());

  dgq_dof_handler.distribute_dofs (fe_dgq);
  dgq_solution.reinit(dgq_dof_handler.n_dofs());

  ConstraintMatrix cm;
  cm.close();
  VectorTools::interpolate (mapping, q_dof_handler, function, q_solution);
  VectorTools::interpolate (mapping, dgq_dof_handler, function, dgq_solution);

  SolutionTransfer<dim,Vector<double>,hp::DoFHandler<dim> > q_soltrans(q_dof_handler);
  SolutionTransfer<dim,Vector<double>,hp::DoFHandler<dim> > dgq_soltrans(dgq_dof_handler);


  // test b): do some coarsening and
  // refinement
  q_soltrans.clear();
  dgq_soltrans.clear();

  counter = 0;
  cell=tria.begin_active();
  endc=tria.end();
  for (; cell!=endc; ++cell, ++counter)
    {
      if (counter > 120)
        cell->set_coarsen_flag();
      else if (Testing::rand() % 3 == 0)
        cell->set_refine_flag();
      else if (Testing::rand() % 3 == 3)
        cell->set_coarsen_flag();
    }

  Vector<double> q_old_solution=q_solution,
                 dgq_old_solution=dgq_solution;
  tria.prepare_coarsening_and_refinement();
  q_soltrans.prepare_for_coarsening_and_refinement(q_old_solution);
  dgq_soltrans.prepare_for_coarsening_and_refinement(dgq_old_solution);
  tria.execute_coarsening_and_refinement();

  counter = 0;
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = q_dof_handler.begin_active(),
    celldg = dgq_dof_handler.begin_active(),
    endc = q_dof_handler.end();
    for (; cell!=endc; ++cell, ++celldg, ++counter)
      {
        if (counter > 20 && counter < 90)
          cell->set_active_fe_index(0);
        else
          cell->set_active_fe_index(Testing::rand()%max_degree);
        if (counter > 20 && counter < 90)
          celldg->set_active_fe_index(0);
        else
          celldg->set_active_fe_index(Testing::rand()%max_degree);
      }
  }

  q_dof_handler.distribute_dofs (fe_q);
  dgq_dof_handler.distribute_dofs (fe_dgq);
  q_solution.reinit(q_dof_handler.n_dofs());
  dgq_solution.reinit(dgq_dof_handler.n_dofs());
  q_soltrans.interpolate(q_old_solution, q_solution);
  dgq_soltrans.interpolate(dgq_old_solution, dgq_solution);

  // check correctness by comparing the values
  // on points of QGauss of order 2.
  MyFunction<dim> func;
  {
    double error = 0;
    const hp::QCollection<dim> quad (QGauss<dim> (2));
    hp::FEValues<dim> hp_fe_val (fe_q, quad, update_values |
                                 update_quadrature_points);
    std::vector<double> vals (quad[0].size());
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = q_dof_handler.begin_active(),
    endc = q_dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        hp_fe_val.reinit (cell, 0);
        const FEValues<dim> &fe_val = hp_fe_val.get_present_fe_values();
        fe_val.get_function_values (q_solution, vals);
        for (unsigned int q=0; q<fe_val.n_quadrature_points; ++q)
          {
            error += std::fabs(func.value(fe_val.quadrature_point(q),0)-
                               vals[q]);
          }
      }
    deallog << "Error in interpolating hp FE_Q: " << error << std::endl;
  }
  {
    double error = 0;
    const hp::QCollection<dim> quad (QGauss<dim> (2));
    hp::FEValues<dim> hp_fe_val (fe_dgq, quad, update_values |
                                 update_quadrature_points);
    std::vector<double> vals (quad[0].size());
    typename hp::DoFHandler<dim>::active_cell_iterator
    celldg = dgq_dof_handler.begin_active(),
    endc =  dgq_dof_handler.end();
    for (; celldg!=endc; ++celldg)
      {
        hp_fe_val.reinit (celldg, 0);
        const FEValues<dim> &fe_val = hp_fe_val.get_present_fe_values();
        fe_val.get_function_values (dgq_solution, vals);
        for (unsigned int q=0; q<fe_val.n_quadrature_points; ++q)
          {
            error += std::fabs(func.value(fe_val.quadrature_point(q),0)-
                               vals[q]);
          }
      }
    deallog << "Error in interpolating hp FE_DGQ: " << error << std::endl;
  }
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



