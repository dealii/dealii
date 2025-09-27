// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that indeed Triangulation::create_triangulation throws an
// exception if we have distorted cells. this test is like the _01
// case except that it reads the mesh from a file

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> coarse_grid(Triangulation<dim>::none, true);

  GridIn<dim>   gi;
  std::ifstream in(dim == 2 ? SOURCE_DIR "/grids/2d" : SOURCE_DIR "/grids/3d");

  gi.attach_triangulation(coarse_grid);


  bool flag = false;
  try
    {
      gi.read_ucd(in);
    }
  catch (ExceptionBase &exc)
    {
      deallog << exc.get_exc_name() << std::endl;
      flag = true;
    }

  Assert(flag == true, ExcInternalError());
}


int
main()
{
  initlog();

  check<2>();
  check<3>();
}
