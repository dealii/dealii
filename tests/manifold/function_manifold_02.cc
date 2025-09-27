// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test the identity Manifold.

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// Helper function
template <int dim, int spacedim>
void
test(unsigned int ref = 1)
{
  deallog << "Testing dim " << dim << ", spacedim " << spacedim << std::endl;

  // Here the only allowed axis is z. In cylinder the default is x.
  std::string push_forward_expression;
  std::string pull_back_expression;

  switch (spacedim)
    {
      case 2:
        push_forward_expression = "x^2; y^2";
        pull_back_expression    = "sqrt(x); sqrt(y)";
        break;
      case 3:
        push_forward_expression = "x^2; y^2; z^2";
        pull_back_expression    = "sqrt(x); sqrt(y); sqrt(z)";
        break;
      default:
        Assert(false, ExcInternalError());
    }

  FunctionManifold<dim, spacedim, spacedim> manifold(push_forward_expression,
                                                     pull_back_expression);

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria, 0, 1);

  for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    {
      cell->set_all_manifold_ids(1);
    }

  tria.set_manifold(1, manifold);
  tria.refine_global(2);

  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();


  test<2, 2>();
  test<3, 3>();

  return 0;
}
