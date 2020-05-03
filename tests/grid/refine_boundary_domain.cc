// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <set>
#include <tuple>

#include "../tests.h"

using namespace dealii;
using namespace std;

int
main()
{
  initlog();

  const unsigned int spacedim = 3;

  // ball as domain triangulation
  Triangulation<spacedim, spacedim> tria_domain;
  GridGenerator::hyper_ball(tria_domain);

  // sphere as boundary triangulation
  Triangulation<spacedim - 1, spacedim> tria_boundary;
  GridGenerator::extract_boundary_mesh(tria_domain, tria_boundary);
  SphericalManifold<spacedim - 1, spacedim> spherical_manifold;
  tria_boundary.set_all_manifold_ids(0);
  tria_boundary.set_manifold(0, spherical_manifold);


  // refine two times
  tria_domain.refine_global(2);
  tria_boundary.refine_global(2);

  // get sets with vertices for boundary cells and domain faces (use float to
  // avoid different ordering in sets caused by round-off)
  set<tuple<const float, const float, const float>> boundary_vertices,
    domain_vertices;
  for (const auto &boundary_cell : tria_boundary.active_cell_iterators())
    {
      for (unsigned int v = 0;
           v < GeometryInfo<spacedim - 1>::vertices_per_cell;
           ++v)
        {
          const Point<spacedim> p = boundary_cell->vertex(v);
          boundary_vertices.insert(make_tuple(p(0), p(1), p(2)));
        }
    }
  for (const auto &domain_cell : tria_domain.active_cell_iterators())
    {
      for (unsigned int f = 0; f < GeometryInfo<spacedim>::faces_per_cell; ++f)
        {
          if (domain_cell->face(f)->at_boundary())
            {
              for (unsigned int v = 0;
                   v < GeometryInfo<spacedim>::vertices_per_face;
                   ++v)
                {
                  const Point<spacedim> p = domain_cell->face(f)->vertex(v);
                  domain_vertices.insert(make_tuple(p(0), p(1), p(2)));
                }
            }
        }
    }

  // compare sets
  if (boundary_vertices != domain_vertices)
    deallog << "ERROR : Meshes do not coincide!" << endl;
  else
    deallog << "OK : Meshes coincide!" << endl;
}
