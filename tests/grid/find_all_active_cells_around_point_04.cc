// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2023 by the deal.II authors
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

// find_all_active_cells_around_point for neighbor cells at different levels.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
print_result(const Mapping<dim, spacedim>       &mapping,
             const Triangulation<dim, spacedim> &tria,
             const Point<dim>                    p)
{
  deallog << "Testing " << dim << "D with point " << p << " tolerance default "
          << std::endl;
  auto c_p = GridTools::find_all_active_cells_around_point(mapping, tria, p);
  for (auto i : c_p)
    deallog << "Cell: " << i.first->id() << " unit point " << i.second
            << std::endl;
  deallog << std::endl;
}

template <int dim, int spacedim>
void
test(const bool flip_refinement = false)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube(
    tria, 2, 0 /*left*/, 1 /*right*/, true /*colorize*/);

  // In every cycle cells at the boundary face x=0 are refined.
  for (unsigned int i = 1; i < 3; ++i)
    {
      for (const auto &cell : tria.active_cell_iterators())
        if (cell->at_boundary(static_cast<int>(flip_refinement)))
          cell->set_refine_flag();

      tria.execute_coarsening_and_refinement();
    }

  MappingQ<dim> mapping(2);

  Point<dim> point_at_edge;

  // boundary face
  point_at_edge[dim - 1] = flip_refinement ? 0.75 : 0.25;
  print_result(mapping, tria, point_at_edge);

  if (dim > 1)
    {
      // interior face
      for (unsigned d = 0; d < dim; ++d)
        point_at_edge[d] = flip_refinement ? 0.75 : 0.25;
      print_result(mapping, tria, point_at_edge);

      // hanging node
      for (unsigned d = 0; d < dim - 1; ++d)
        point_at_edge[d] = 0.5;
      print_result(mapping, tria, point_at_edge);
    }
  else
    {
      // interior face
      point_at_edge[0] = 0.5;
      print_result(mapping, tria, point_at_edge);
    }

  // debug
  if (false)
    {
      std::ofstream output_file("grid_" + std::to_string(dim) + ".vtu");
      GridOut().write_vtu(tria, output_file);
    }
};


int
main()
{
  initlog();
  deallog << std::setprecision(10);
  test<1, 1>();
  test<2, 2>();
  test<3, 3>();
  test<1, 1>(true);
  test<2, 2>(true);
  test<3, 3>(true);
  return 0;
}
