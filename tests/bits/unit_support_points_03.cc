// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test FiniteElement::unit_support_points() and ::unit_face_support_points()
// for different shapes.

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_wedge_p.h>

#include "../tests.h"


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
