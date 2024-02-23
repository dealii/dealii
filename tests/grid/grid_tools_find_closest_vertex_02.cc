// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check find_closest_vertex_of_cell when a mapping is present

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test(const Point<spacedim> &p, double displacement)
{
  deallog << "dim: " << dim << ", spacedim: " << spacedim
          << " with displacement " << displacement << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global();

  // Vector FE
  FESystem<dim, spacedim>   fe{FE_Q<dim, spacedim>(1), spacedim};
  DoFHandler<dim, spacedim> dh(tria);
  dh.distribute_dofs(fe);

  // Displacement vector.
  Vector<double> disp(dh.n_dofs());
  for (unsigned int i = 0; i < dh.n_dofs(); ++i)
    disp[i] = displacement;

  // mapping
  MappingQEulerian<dim, Vector<double>, spacedim> mapping(1, dh, disp);

  for (const auto cell : tria.active_cell_iterators())
    {
      const auto i =
        GridTools::find_closest_vertex_of_cell<dim, spacedim>(cell, p, mapping);
      const auto &v = mapping.get_vertices(cell);

      deallog << "Closest vertex of cell " << cell << " to " << p
              << " is cell->vertex(" << i << ") : " << v[i] << std::endl;
    }
};


int
main()
{
  initlog();
  test<2, 2>(Point<2>(.45, .45), 0);
  test<2, 2>(Point<2>(.45, .45), 0.5);

  test<2, 3>(Point<3>(.9, .9, .1), 0.0);
  test<2, 3>(Point<3>(.9, .9, .1), 1.0);
  return 0;
}
