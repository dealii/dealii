// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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


// find_all_active_cells_around_point for a simple cube mesh

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/mapping_q.h>

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
test(unsigned int n_ref)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_ref);
  MappingQ<dim, spacedim> mapping(3);

  Point<dim> p;
  {
    for (unsigned int d = 0; d < dim; ++d)
      p[d] = 0.5;
    print_result(mapping, tria, p, 1e-8);
  }
  {
    for (unsigned int d = 0; d < dim; ++d)
      p[d] = 0.5 + 1e-6;
    print_result(mapping, tria, p, 1e-12);
    print_result(mapping, tria, p, 1e-4);
  }
  for (unsigned int c = 0; c < dim; ++c)
    for (unsigned int e = 0; e < dim; ++e)
      {
        for (unsigned int d = 0; d < dim; ++d)
          p[d] = 0.5;
        p[c] = 0.45;
        p[e] = 0.45;
        print_result(mapping, tria, p, 1e-12);
      }
}



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  test<1, 1>(4);
  test<2, 2>(4);
  test<2, 2>(5);
  test<3, 3>(3);
}
