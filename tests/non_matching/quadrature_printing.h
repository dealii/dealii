// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2021 by the deal.II authors
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

#ifndef dealii_quadrature_printing_h_
#define dealii_quadrature_printing_h_

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/non_matching/immersed_surface_quadrature.h>

using namespace dealii;

/*
 * Print the incoming quadrature to deallog as comma separated values:
 * point[0], ..., point[dim-1], weight
 */
template <int dim>
void
print_quadrature(const Quadrature<dim> &quadrature)
{
  for (unsigned int i = 0; i < quadrature.size(); ++i)
    {
      const Point<dim> &point = quadrature.point(i);
      for (int d = 0; d < dim; ++d)
        deallog << point[d] << ", ";

      deallog << quadrature.weight(i) << std::endl;
    }
}



/*
 * Print the incoming surface quadrature to deallog as comma separated values:
 * p[0], ..., p[dim-1], weight, normal[0], ..., normal[spacedim-1]
 */
template <int dim, int spacedim>
void
print_surface_quadrature(
  const NonMatching::ImmersedSurfaceQuadrature<dim, spacedim> &quadrature)
{
  for (unsigned int i = 0; i < quadrature.size(); ++i)
    {
      const Point<dim> &point = quadrature.point(i);
      for (int d = 0; d < dim; ++d)
        deallog << point[d] << ", ";

      deallog << quadrature.weight(i);

      const Tensor<1, spacedim> &normal = quadrature.normal_vector(i);
      for (int d = 0; d < spacedim; ++d)
        deallog << ", " << normal[d];
      deallog << std::endl;
    }
}


#endif
