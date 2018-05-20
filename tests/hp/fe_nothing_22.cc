// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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

// Test FE_Nothing::operator==()

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/numerics/vector_tools.h>

template <int dim>
void
test()
{
  deallog << std::boolalpha;
  deallog << (FE_Nothing<dim>(1) == FE_Nothing<dim>(1, false)) << std::endl;
  deallog << (FE_Nothing<dim>(1) == FE_Nothing<dim>(2)) << std::endl;
  deallog << (FE_Nothing<dim>(2, true) == FE_Nothing<dim>(2, false))
          << std::endl;
  deallog << (FE_Nothing<dim>(1, true) == FE_Nothing<dim>(2, true))
          << std::endl;
}

int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;
}
