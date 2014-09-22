// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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

// document a bug in hp::DoFHandler and copy_triangulation

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/hp/dof_handler.h>

#include <fstream>
#include <cmath>

using namespace dealii;



template <int dim>
void test1()
{
  Triangulation<dim> tr1;
  hp::DoFHandler<dim> dofh(tr1);

  Triangulation<dim> tr2;
  GridGenerator::hyper_ball(tr2);
  tr2.refine_global(2);
  tr1.copy_triangulation(tr2);

  dofh.begin_active()->set_active_fe_index(0);
  
  
}

int main()
{
  initlog();

  test1<2>();
  return 0;
}
