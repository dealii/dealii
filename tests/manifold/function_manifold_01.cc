// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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

// Test the identity Manifold.

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

// Helper function
template <int dim, int spacedim>
void
test(unsigned int ref=1)
{
  deallog << "Testing dim " << dim
          << ", spacedim " << spacedim << std::endl;

  // Here the only allowed axis is z. In cylinder the default is x.
  std::string push_forward_expression;
  std::string pull_back_expression;

  switch (spacedim)
    {
    case 2:
      push_forward_expression = "x; y";
      pull_back_expression = "x; y";
      break;
    case 3:
      push_forward_expression = "x; y; z";
      pull_back_expression = "x; y; z";
      break;
    default:
      Assert(false, ExcInternalError());
    }

  FunctionManifold<dim,spacedim,spacedim> manifold(push_forward_expression,
                                                   pull_back_expression);

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube (tria, 0, 1);

  for (typename Triangulation<dim,spacedim>::active_cell_iterator cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      cell->set_all_manifold_ids(1);
    }

  tria.set_manifold(1, manifold);
  tria.refine_global(2);

  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
}

int
main ()
{
  initlog();


  test<2,2>();
  test<3,3>();

  return 0;
}
