// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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
