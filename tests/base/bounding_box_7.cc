// ------------------------------------------------------------------------
//
// Copyright (C) 2020 - 2024 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test the following functions of the BoundingBox class
// unit_to_real
// real_to_unit

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/std_cxx20/iota_view.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>

#include "../tests.h"


using iota = std_cxx20::ranges::iota_view<unsigned int, unsigned int>;

template <int dim>
void
test()
{
  const unsigned int n_boxes  = 3;
  const unsigned int n_points = 3;

  const auto unit_box = create_unit_bounding_box<dim>();

  for (auto i : iota(0, n_boxes))
    {
      const auto box = random_box<dim>();
      // Check that all vertices get transformed correctly
      for (auto v : iota(0, GeometryInfo<dim>::vertices_per_cell))
        {
          const auto vertex      = unit_box.vertex(v);
          const auto real_vertex = box.unit_to_real(vertex);
          if (real_vertex.distance(box.vertex(v)) > 1e-10)
            deallog << "Error: " << vertex << " -> " << real_vertex
                    << " != " << box.vertex(v) << std::endl;

          const auto unit_vertex = box.real_to_unit(real_vertex);
          if (unit_vertex.distance(vertex) > 1e-10)
            deallog << "Error: " << real_vertex << " -> " << unit_vertex
                    << " != " << vertex << std::endl;
        }

      // Now check random points
      for (auto j : iota(0, n_points))
        {
          const auto unit_point   = random_point<dim>();
          const auto real_point   = box.unit_to_real(unit_point);
          const auto real_to_unit = box.real_to_unit(real_point);

          if (real_to_unit.distance(unit_point) > 1e-10)
            deallog << "Error: " << unit_point << " -> " << real_point << " -> "
                    << real_to_unit << " != " << unit_point << std::endl;
        }
    }
  deallog << "Spacedim " << dim << " OK" << std::endl;
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
