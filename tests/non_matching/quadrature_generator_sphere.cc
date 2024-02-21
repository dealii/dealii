// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include <cmath>

#include "../tests.h"


/**
 * Compute the volume and surface area of a ball/sphere by setting up a
 * level set function immersed in a background mesh, generating the
 * quadrature rules, and summing the weights.
 */
template <int dim>
void
calculate_volume_and_surface_area()
{
  // Set up a background mesh
  Triangulation<dim> triangulation;
  const int          n_subdivisions = 12;
  const double       gridsize       = 2.07;
  GridGenerator::subdivided_hyper_cube(triangulation,
                                       n_subdivisions,
                                       -gridsize / 2,
                                       gridsize / 2);


  // Description of the immersed domain.
  const Functions::SignedDistance::Sphere<dim> level_set;

  // Create a quadrature generator.
  const hp::QCollection<1>              q_collection1D(QGauss<1>(2));
  NonMatching::QuadratureGenerator<dim> quadrature_generator(q_collection1D);

  // Go over all cells and compute the volume and surface area.
  double surface_area = 0, volume = 0;
  for (const auto cell : triangulation.active_cell_iterators())
    {
      // Create a box corresponding to the cell.
      std::pair<Point<dim>, Point<dim>> lower_upper_corner;
      lower_upper_corner.first = cell->vertex(0);
      lower_upper_corner.second =
        cell->vertex(GeometryInfo<dim>::vertices_per_cell - 1);
      const BoundingBox<dim> box(lower_upper_corner);

      // Generate immersed quadrature rules.
      quadrature_generator.generate(level_set, box);

      // Get the quadrature rules.
      const Quadrature<dim> &inside_quadrature =
        quadrature_generator.get_inside_quadrature();
      const NonMatching::ImmersedSurfaceQuadrature<dim> &surface_quadrature =
        quadrature_generator.get_surface_quadrature();

      // Sum the weights to get the area/volume of the sphere.
      for (unsigned int i = 0; i < inside_quadrature.size(); ++i)
        volume += inside_quadrature.weight(i);

      // Sum the weights to get the circumference/surface area of the sphere.
      for (unsigned int i = 0; i < surface_quadrature.size(); ++i)
        surface_area += surface_quadrature.weight(i);
    }

  deallog << "dim = " << dim << std::endl;
  deallog << (2 == dim ? "area = " : "volume = ");
  deallog << volume / M_PI << " * pi" << std::endl;
  deallog << (2 == dim ? "circumference = " : "surface area = ");
  deallog << surface_area / M_PI << " * pi" << std::endl;
  deallog << std::endl;
}


int
main()
{
  initlog();
  calculate_volume_and_surface_area<2>();
  calculate_volume_and_surface_area<3>();
}
