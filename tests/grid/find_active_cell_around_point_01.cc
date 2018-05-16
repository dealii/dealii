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


/*
 * Given the number of refinements and the number of random points
 * it benchmarks the time needed to run the function FCT
 * which can be point_locator_D2 (or point_locator when it shall be written)
 */
#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>


template <int dim, int spacedim>
void
test(unsigned int n_ref, unsigned int n_points)
{
  deallog << "Testing " << dim << ", " << spacedim << std::endl;

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_ref);
  DoFHandler<dim,spacedim> dof_handler(tria);

  std::vector<Point<spacedim>> points;

  deallog << "Points in study: " << n_points << std::endl;
  for (size_t i=0; i<n_points; ++i)
    points.push_back(random_point<spacedim>());

  auto v_to_c = GridTools::vertex_to_cell_map (tria);
  auto v_to_c_d = GridTools::vertex_to_cell_centers_directions (tria, v_to_c);

  auto &mapping = StaticMappingQ1<dim,spacedim>::mapping;
  auto cell = tria.begin_active();
  for (auto &p : points)
    {
      auto c_and_p = GridTools::find_active_cell_around_point(mapping, tria, p,
                                                              v_to_c, v_to_c_d);
      auto p2 = mapping.transform_unit_to_real_cell(c_and_p.first, c_and_p.second);
      if (p2.distance(p)>1e-10)
        deallog << "NOT OK!" << p << ", " << p2
                << ", " << c_and_p.first << std::endl;
      cell = c_and_p.first;
    }
  deallog << "OK" << std::endl;
}

int
main ()
{
  initlog();

  test<1,1> (4,10);
  test<2,2> (3,20);
  test<3,3> (2,30);
}
