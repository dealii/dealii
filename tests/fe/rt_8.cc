//----------------------------  rt_8.cc  ---------------------------
//    rt_8.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2003, 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rt_8.cc  ---------------------------

// build a mass matrix for the RT element and try to invert it. we had trouble
// with this at one time

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/vector_memory.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 5


std::ofstream logfile ("rt_8/output");

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


  for (unsigned int q_point=0; q_point<q.n_quadrature_points; ++q_point)
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
        for (unsigned int d=0; d<dim; ++d)
          mass_matrix(i,j)
            += (fe.shape_value_component(i,q_point,d) *
                fe.shape_value_component(j,q_point,d) *
                fe.JxW(q_point));
  mass_matrix.print_formatted (logfile, 3, false, 0, " ", 1);

  SolverControl           solver_control (dofs_per_cell,
                                          1e-8);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);

  Vector<double> tmp1(dofs_per_cell), tmp2(dofs_per_cell);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    tmp1(i) = 1.*rand()/RAND_MAX;
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



