// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// VectorTools::interpolate_boundary_values still had bugs in 1d after
// switching to a scheme where we can assign boundary indicators also in 1d

#include "../tests.h"
#include <deal.II/base/function_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iomanip>
#include <vector>



template <int dim>
void test ()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim>   triangulation;
  FE_Q<dim>            fe(1);
  DoFHandler<dim>      dof_handler(triangulation);

  GridGenerator::hyper_cube (triangulation, 0, 1);
  triangulation.begin_active()->face(0)->set_boundary_indicator(10);
  triangulation.begin_active()->face(1)->set_boundary_indicator(20);
  triangulation.refine_global (1);

  dof_handler.distribute_dofs (fe);

  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            10,
                                            Functions::SquareFunction<dim>(),
                                            boundary_values);
  VectorTools::interpolate_boundary_values (dof_handler,
                                            20,
                                            Functions::SquareFunction<dim>(),
                                            boundary_values);
  deallog << boundary_values.size() << std::endl;
  for (std::map<types::global_dof_index,double>::const_iterator
       p = boundary_values.begin();
       p != boundary_values.end(); ++p)
    deallog << p->first << ' ' << p->second << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
}

