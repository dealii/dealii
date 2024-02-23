// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// tests which used to fail because influence of vertex swapping (to get proper
// normals) on numbering of children was not taken into account in
// extract_boundary_mesh

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>

#include "../tests.h"


template <int spacedim>
void
test(std::ostream &out)
{
  out << "spacedim=" << spacedim << std::endl << std::endl;

  // triangulation of domain (unit cube refined once globally)
  Triangulation<spacedim> triaDomain;
  GridGenerator::hyper_cube(triaDomain, -1, 1);

  // refinement
  unsigned int refined_cell;
  if (spacedim == 2)
    refined_cell = 3;
  else if (spacedim == 3)
    refined_cell = 3;
  else
    return;

  auto cell = triaDomain.begin(0);
  triaDomain.refine_global(1);
  cell->child(refined_cell)->set_refine_flag();
  triaDomain.execute_coarsening_and_refinement();
  cell->child(refined_cell)->child(refined_cell)->set_refine_flag();
  triaDomain.execute_coarsening_and_refinement();

  // extract triangulation of boundary
  Triangulation<spacedim - 1, spacedim> tria_boundary;
  auto                                  map_boundary_to_domain =
    GridGenerator::extract_boundary_mesh(triaDomain, tria_boundary);

  // output of vertex positions
  out << "Vertex positions, center difference:" << std::endl;
  for (auto cell_pair = map_boundary_to_domain.begin();
       cell_pair != map_boundary_to_domain.end();
       ++cell_pair)
    {
      for (unsigned int vertex = 0;
           vertex < GeometryInfo<spacedim>::vertices_per_face;
           vertex++)
        out << cell_pair->first->vertex(vertex) << std::endl;
      out << "-----------------------" << std::endl;
      for (unsigned int vertex = 0;
           vertex < GeometryInfo<spacedim>::vertices_per_face;
           vertex++)
        out << cell_pair->second->vertex(vertex) << std::endl;
      out << "-----------------------" << std::endl;
      out << cell_pair->first->center() - cell_pair->second->center()
          << std::endl
          << std::endl;
    }

  // output of boundary triangulation
  GridOut go;
  out << std::endl << "Gnuplot:" << std::endl;
  go.write_gnuplot(tria_boundary, out);
}

int
main()
{
  initlog();
  test<2>(deallog.get_file_stream());
  test<3>(deallog.get_file_stream());
}
