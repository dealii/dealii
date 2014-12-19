// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


// build a mass matrix for the RT element and try to invert it. we had trouble
// with this at one time

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 5


std::ofstream logfile ("output");

template<int dim>
void
test (const unsigned int degree)
{
  FE_RaviartThomas<dim> fe_rt(degree);
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

  DoFHandler<dim> dof(tr);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(fe_rt);

  QTrapez<1> q_trapez;
  const unsigned int div=4;
  QIterated<dim> q(q_trapez, div);
  FEValues<dim> fe(fe_rt, q, update_values|update_JxW_values);
  fe.reinit(c);

  const unsigned int dofs_per_cell = fe_rt.dofs_per_cell;
  FullMatrix<double> mass_matrix (dofs_per_cell, dofs_per_cell);

  Assert (fe.get_fe().n_components() == dim, ExcInternalError());


  for (unsigned int q_point=0; q_point<q.size(); ++q_point)
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
        for (unsigned int d=0; d<dim; ++d)
          mass_matrix(i,j)
          += (fe.shape_value_component(i,q_point,d) *
              fe.shape_value_component(j,q_point,d) *
              fe.JxW(q_point));
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    for (unsigned int j=0; j<dofs_per_cell; ++j)
      if (std::fabs(mass_matrix(i,j)) < 1e-14)
        mass_matrix(i,j) = 0;
  mass_matrix.print_formatted (logfile, 3, false, 0, " ", 1);

  SolverControl           solver_control (dofs_per_cell,
                                          1e-8);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);

  Vector<double> tmp1(dofs_per_cell), tmp2(dofs_per_cell);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    tmp1(i) = 1.*Testing::rand()/RAND_MAX;
  cg.solve (mass_matrix, tmp2, tmp1, PreconditionIdentity());

  deallog << "Degree=" << degree
          << ": " << solver_control.last_step()
          << " iterations to obtain convergence."
          << std::endl;
}


int
main()
{
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int i=0; i<4; ++i)
    test<2>(i);

  return 0;
}



