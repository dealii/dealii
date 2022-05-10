// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

// Test FiniteElement::unit_support_points() and ::unit_face_support_points()
// for different shapes.

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"

using namespace dealii;

template <int dim, int spacedim>
void
test(const FiniteElement<dim, spacedim> &fe)
{
  deallog << fe.get_name() << std::endl;

  const auto &unit_support_points = fe.get_unit_support_points();

  AssertDimension(unit_support_points.size(), fe.n_dofs_per_cell());

  for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
    for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
      {
        const auto value = fe.shape_value(i, unit_support_points[j]);
        Assert((i == j ? (std::abs(value - 1.0) < 1e-8) :
                         (std::abs(value) < 1e-8)),
               ExcInternalError());
      }

  FE_Q<dim - 1>        fe_q(fe.degree);
  FE_SimplexP<dim - 1> fe_p(fe.degree);

  for (const auto f : fe.reference_cell().face_indices())
    {
      const auto &unit_face_support_points = fe.get_unit_face_support_points(f);

      AssertDimension(unit_face_support_points.size(), fe.n_dofs_per_face(f));

      const auto &fe_face =
        fe.reference_cell().face_reference_cell(f).is_hyper_cube() ?
          static_cast<FiniteElement<dim - 1> &>(fe_q) :
          static_cast<FiniteElement<dim - 1> &>(fe_p);

      for (unsigned int i = 0; i < fe.n_dofs_per_face(f); ++i)
        for (unsigned int j = 0; j < fe.n_dofs_per_face(f); ++j)
          {
            const auto value =
              fe_face.shape_value(i, unit_face_support_points[j]);
            Assert((i == j ? (std::abs(value - 1.0) < 1e-8) :
                             (std::abs(value) < 1e-8)),
                   ExcInternalError());
          }
    }



  deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  test(FE_Q<2>(1));
  test(FE_Q<2>(2));
  test(FE_Q<3>(1));
  test(FE_Q<3>(2));

  test(FE_SimplexP<2>(1));
  test(FE_SimplexP<2>(2));
  test(FE_SimplexP<3>(1));
  test(FE_SimplexP<3>(2));

  test(FE_PyramidP<3>(1));

  test(FE_WedgeP<3>(1));
  test(FE_WedgeP<3>(2));
}
