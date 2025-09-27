// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check extract used_vertices and find_closest_vertex when a mapping
// is present

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
test(const Point<spacedim> &p)
{
  deallog << "dim: " << dim << ", spacedim: " << spacedim << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  // Vector FE
  FESystem<dim, spacedim>   fe{FE_Q<dim, spacedim>(1), spacedim};
  DoFHandler<dim, spacedim> dh(tria);
  dh.distribute_dofs(fe);

  // Displacement vector.
  Vector<double> disp(dh.n_dofs());
  for (unsigned int i = 0; i < dh.n_dofs(); ++i)
    disp[i] = .5;

  // mapping
  MappingQEulerian<dim, Vector<double>, spacedim> mapping(1, dh, disp);

  auto m = GridTools::extract_used_vertices(tria, mapping);

  std::vector<bool> selected_vertices(tria.n_vertices(), false);
  selected_vertices[0] = true;
  selected_vertices[5] = true;

  for (auto &e : m)
    deallog << "Vertex: " << e.first << ": " << e.second << std::endl;

  auto i = GridTools::find_closest_vertex(mapping, tria, p, selected_vertices);
  deallog << "Closest vertex to " << p << ", v[" << i << "] :" << m[i]
          << std::endl;
};


int
main()
{
  initlog();
  test<2, 2>(Point<2>(.2, .2));
  test<2, 2>(Point<2>(.6, .9));
  return 0;
}
