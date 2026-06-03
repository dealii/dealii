// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check GridTools::affine_cell_approximation for non-hypercube reference cells.

#include <deal.II/base/array_view.h>
#include <deal.II/base/derivative_form.h>

#include <deal.II/grid/grid_tools_geometry.h>
#include <deal.II/grid/reference_cell.h>

#include "../tests.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>


template <int dim, int spacedim>
DerivativeForm<1, dim, spacedim>
make_linear_part()
{
  DerivativeForm<1, dim, spacedim> A;
  for (unsigned int d = 0; d < spacedim; ++d)
    for (unsigned int e = 0; e < dim; ++e)
      A[d][e] = (d == e ? 1.0 + 0.2 * d : 0.0) +
                0.05 * (d + 1) * (e + 2);

  return A;
}



template <int spacedim>
Tensor<1, spacedim>
make_offset()
{
  Tensor<1, spacedim> b;
  for (unsigned int d = 0; d < spacedim; ++d)
    b[d] = -0.3 + 0.2 * d;

  return b;
}



template <int dim, int spacedim>
Point<spacedim>
transform(const DerivativeForm<1, dim, spacedim> &A,
          const Tensor<1, spacedim>              &b,
          const Point<dim>                       &unit_point)
{
  Point<spacedim> real_point;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      real_point[d] = b[d];
      for (unsigned int e = 0; e < dim; ++e)
        real_point[d] += A[d][e] * unit_point[e];
    }

  return real_point;
}



template <int dim, int spacedim>
void
test(const ReferenceCell<dim> &reference_cell, const std::string &name)
{
  const auto A = make_linear_part<dim, spacedim>();
  const auto b = make_offset<spacedim>();

  std::vector<Point<spacedim>> vertices(reference_cell.n_vertices());
  for (unsigned int v = 0; v < reference_cell.n_vertices(); ++v)
    vertices[v] = transform(A, b, reference_cell.vertex(v));

  const auto A_b =
    GridTools::affine_cell_approximation<dim, spacedim>(
      reference_cell, make_array_view(vertices));

  double error = 0.0;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      error = std::max(error, std::fabs(A_b.second[d] - b[d]));
      for (unsigned int e = 0; e < dim; ++e)
        error = std::max(error, std::fabs(A_b.first[d][e] - A[d][e]));
    }

  AssertThrow(error < 1e-12, ExcMessage("Unexpected affine fit error."));
  deallog << name << " OK" << std::endl;
}



int
main()
{
  initlog();

  test<2, 2>(ReferenceCells::Triangle, "Triangle");
  test<3, 3>(ReferenceCells::Tetrahedron, "Tetrahedron");
  test<3, 3>(ReferenceCells::Wedge, "Wedge");
  test<3, 3>(ReferenceCells::Pyramid, "Pyramid");
}
