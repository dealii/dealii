// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// check that assigning one invalid iterator to another works. this
// test is for ::DoFHandler

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
void
check()
{
  typename DoFHandler<dim>::active_cell_iterator invalid_1;

  // try copy constructor
  typename DoFHandler<dim>::active_cell_iterator invalid_2 = invalid_1;

  // now also try copy operator
  invalid_1 = invalid_2;

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
}
