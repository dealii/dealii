// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that MappingQCache based on MappingQ1 works correctly in case some
// deformation is applied and CellSimilarity would get activated for the
// undeformed geometry, but not in the deformed state

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
do_test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  tria.refine_global(4 - dim);

  MappingQ<dim>      mapping(1);
  MappingQCache<dim> mapping_cache(1);
  const auto transform = [](const typename Triangulation<dim>::cell_iterator &,
                            const Point<dim> &in) {
    Point<dim> out;
    // move points more densely near zero
    for (unsigned int d = 0; d < dim; ++d)
      out[d] = std::pow(in[d], 1.2 + 0.2 * d);
    return out;
  };
  mapping_cache.initialize(mapping, tria, transform, false);

  const Vector<double> aspect_ratios =
    GridTools::compute_aspect_ratio_of_cells(mapping, tria, QGauss<dim>(2));
  const Vector<double> aspect_ratios_deformed =
    GridTools::compute_aspect_ratio_of_cells(mapping_cache,
                                             tria,
                                             QGauss<dim>(2));

  deallog.push(std::to_string(dim) + "d");
  deallog << "Aspect ratios undeformed mesh: " << std::endl;
  aspect_ratios.print(deallog.get_file_stream());
  deallog << std::endl;
  deallog << "Aspect ratios deformed mesh: " << std::endl;
  aspect_ratios_deformed.print(deallog.get_file_stream());
  deallog << std::endl;
  deallog.pop();
  deallog << std::endl;
}


int
main()
{
  initlog();
  MultithreadInfo::set_thread_limit(1);
  do_test<1>();
  do_test<2>();
  do_test<3>();
}
