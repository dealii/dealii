// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <set>
#include <tuple>

#include "../tests.h"


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

  // create a set of boundary and domain vertices (that have to coincide).
  // In order to compare the two sets we order them in lexicographic order
  // and identify coordinates as being identical up to a certain tolerance,
  // i.e., is_less_tol only returns true if one (or more) coordinates of
  // two points have a difference of at least tol.
  //
  // Warning: Ensure that is_less_tol is a "weak strict ordering", i.e.
  // ensure that is_less_tol(p, p) == false for any p.

  const double tol = 1.0e-8;

  auto is_less_tol = [tol](const Point<spacedim> &p_1,
                           const Point<spacedim> &p_2) {
    if (p_1[0] < p_2[0] - tol)
      return true;

    if (p_1[0] <= p_2[0] + tol && p_1[1] < p_2[1] - tol)
      return true;

    if (p_1[0] <= p_2[0] + tol && p_1[1] <= p_2[1] + tol &&
        p_1[2] < p_2[2] - tol)
      return true;

    return false;
  };

  std::set<Point<spacedim>, decltype(is_less_tol)> boundary_vertices(
    is_less_tol);

  for (const auto &boundary_cell : tria_boundary.active_cell_iterators())
    for (unsigned int v = 0; v < GeometryInfo<spacedim - 1>::vertices_per_cell;
         ++v)
      boundary_vertices.insert(boundary_cell->vertex(v));

  std::set<Point<spacedim>, decltype(is_less_tol)> domain_vertices(is_less_tol);

  for (const auto &domain_cell : tria_domain.active_cell_iterators())
    for (unsigned int f = 0; f < GeometryInfo<spacedim>::faces_per_cell; ++f)
      if (domain_cell->face(f)->at_boundary())
        for (unsigned int v = 0; v < GeometryInfo<spacedim>::vertices_per_face;
             ++v)
          {
            domain_vertices.insert(domain_cell->face(f)->vertex(v));
          }

  // Create the symmetric set difference:

  std::vector<Point<spacedim>> symmetric_difference;
  std::set_symmetric_difference(boundary_vertices.begin(),
                                boundary_vertices.end(),
                                domain_vertices.begin(),
                                domain_vertices.end(),
                                std::back_inserter(symmetric_difference),
                                is_less_tol);

  if (!symmetric_difference.empty())
    {
      deallog << "Found nonmatching vertices" << std::endl;
      for (auto it : symmetric_difference)
        deallog << it << std::endl;
    }
  else
    {
      deallog << "OK : Meshes coincide!" << std::endl;
    }
}
