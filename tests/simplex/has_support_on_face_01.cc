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


// Verify FE_SimplexPoly::has_support_on_face() by comparing its answer with
// the values of the shape functions on the faces of the reference cell.

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>

#include "../tests.h"


template <int dim>
bool
vanishes_on_face(const FiniteElement<dim> &fe,
                 const unsigned int        shape_index,
                 const unsigned int        face_index)
{
  const ReferenceCell<dim> cell       = fe.reference_cell();
  const unsigned int       n_vertices = dim == 1 ? 1 : dim;
  const unsigned int       n_steps    = 8;

  std::vector<Point<dim>> vertices;
  for (unsigned int v = 0; v < n_vertices; ++v)
    vertices.push_back(cell.face_vertex_location(face_index, v));

  // Sample the face with convex combinations of its vertices:
  for (unsigned int a = 0; a <= n_steps; ++a)
    for (unsigned int b = 0; b <= (dim == 3 ? n_steps - a : 0); ++b)
      {
        const double w_1 = (dim == 1 ? 0. : 1. * a / n_steps);
        const double w_2 = (dim == 3 ? 1. * b / n_steps : 0.);

        Point<dim> p = (1. - w_1 - w_2) * vertices[0];
        if (dim >= 2)
          p += w_1 * vertices[1];
        if (dim == 3)
          p += w_2 * vertices[n_vertices - 1];

        if (std::abs(fe.shape_value(shape_index, p)) > 1e-10)
          return false;
      }

  return true;
}


template <int dim>
void
check(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  for (unsigned int f = 0; f < fe.reference_cell().n_faces(); ++f)
    {
      deallog << "face " << f << ':';
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        {
          const bool support = fe.has_support_on_face(i, f);
          AssertThrow(support != vanishes_on_face(fe, i, f),
                      ExcInternalError());
          deallog << ' ' << support;
        }
      deallog << std::endl;
    }
}


template <int dim>
void
test()
{
  for (unsigned int degree = 1; degree <= (dim == 1 ? 2 : 3); ++degree)
    check(FE_SimplexP<dim>(degree));

  for (unsigned int degree = 1; degree <= (dim == 1 ? 2 : 3); ++degree)
    check(FE_SimplexDGP<dim>(degree));

  if (dim > 1)
    for (unsigned int degree = 1; degree <= 2; ++degree)
      check(FE_SimplexP_Bubbles<dim>(degree));
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
