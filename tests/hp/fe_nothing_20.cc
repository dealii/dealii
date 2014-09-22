// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// interpolate() can not deal with FE_Nothing in an hp setting


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/function.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/dof_handler.h>


#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim>       triangulation;
  GridGenerator :: hyper_cube (triangulation, -0.5, 0.5);
  triangulation.refine_global(4);

  hp::FECollection<dim>    fe_collection;
  
  fe_collection.push_back (FESystem<dim>(FE_Q<dim>(2),dim,
					 FE_Q<dim>(1),1,
					 FE_Nothing<dim>(), dim));

  fe_collection.push_back (FESystem<dim>(FE_Nothing<dim>(dim+1), 1,
                                         FE_Q<dim>(2), dim));

  hp::DoFHandler<dim>      dof_handler (triangulation);



  dof_handler.distribute_dofs (fe_collection);  

  deallog << "   Number of active cells:       "
          << triangulation.n_active_cells()
          << std::endl
          << "   Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;

  
  Vector<double> solution(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler,
                             ZeroFunction<dim>(2*dim+1),
                             solution);

  deallog << "l2_norm = " << solution.l2_norm() << std::endl;
}



int main ()
{
  initlog();

  //test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
