// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// find_all_active_cells_around_point around a cylinder

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include "../tests.h"



template <int dim, int spacedim>
void
print_result(const Mapping<dim, spacedim> &      mapping,
             const Triangulation<dim, spacedim> &tria,
             const Point<dim>                    p,
             const double                        tolerance)
{
  deallog << "Testing " << dim << "D with point " << p << " tolerance "
          << tolerance << std::endl;
  auto c_p =
    GridTools::find_all_active_cells_around_point(mapping, tria, p, tolerance);
  for (auto i : c_p)
    deallog << "Cell: " << i.first->id() << " unit point " << i.second
            << std::endl;
  deallog << std::endl;
}



template <int dim, int spacedim>
void
print_result(const Mapping<dim, spacedim> &      mapping,
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
test(unsigned int n_ref)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::channel_with_cylinder(tria, 0.03, 2, 2);
  tria.refine_global(n_ref);
  MappingQGeneric<dim, spacedim> mapping(3);

  Point<dim> p1;
  p1[0] = 0.28;
  p1[1] = 0.2;
  if (dim == 3)
    p1[2] = 0.205;
  print_result(mapping, tria, p1, 1e-3);
  print_result(mapping, tria, p1, 1e-4);
  print_result(mapping, tria, p1, 1e-5);
  print_result(mapping, tria, p1, 1e-6);
  print_result(mapping, tria, p1, 1e-7);
  print_result(mapping, tria, p1, 1e-10);
  print_result(mapping, tria, p1);
  Point<dim> p2 = p1;
  // this point lies on exactly 2 cells in 2D and 4 in 3D, up to a
  p2[0] = 0.0682926829268293;
  print_result(mapping, tria, p2, 1e-3);
  print_result(mapping, tria, p2, 1e-4);
  print_result(mapping, tria, p2, 1e-5);
  print_result(mapping, tria, p2, 1e-6);
  print_result(mapping, tria, p2, 1e-7);
  print_result(mapping, tria, p2, 1e-10);
  print_result(mapping, tria, p2);
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<2, 2>(1);
  test<2, 2>(2);
  test<3, 3>(0);
  test<3, 3>(1);
}
