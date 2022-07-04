// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Short test to validate compute_bounding_box()
//
// like _01, but assign the result to a BoundingBox object

#include <deal.II/base/bounding_box.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
bool
pred_mat_id(const typename Triangulation<dim>::active_cell_iterator &cell)
{
  return cell->material_id() == 2;
}

template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  using cell_iterator = typename Triangulation<dim>::active_cell_iterator;

  // Mark a small block at the corner of the hypercube
  cell_iterator cell = tria.begin_active(), endc = tria.end();
  for (; cell != endc; ++cell)
    {
      bool mark = true;
      for (unsigned int d = 0; d < dim; ++d)
        if (cell->center()[d] > 0.33)
          {
            mark = false;
            break;
          }

      if (mark == true)
        cell->set_material_id(2);
      else
        cell->set_material_id(1);
    }

  std::function<bool(const cell_iterator &)> predicate = pred_mat_id<dim>;

  // Find bounding box that surrounds cells with material id 2
  BoundingBox<dim> bounding_box =
    GridTools::compute_bounding_box(tria, predicate); // General predicate

  deallog << bounding_box.get_boundary_points().first << ' '
          << bounding_box.get_boundary_points().second << std::endl;
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
