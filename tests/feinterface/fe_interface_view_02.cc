// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test FEInterfaceViews and extractor classes. these tests use a primitive
// finite element and vector extractors. To simplify the output, we skip writing
// results of high-order tensors and zeros.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
test(const Triangulation<dim> &tr, const FiniteElement<dim> &fe)
{
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  deallog << "FE=" << fe.get_name() << std::endl;

  const QGauss<dim - 1>  face_quadrature(2);
  FEInterfaceValues<dim> fe_iv(fe,
                               face_quadrature,
                               update_values | update_gradients |
                                 update_quadrature_points | update_hessians |
                                 update_3rd_derivatives);
  auto                   cell = dof.begin_active();


  for (const auto f : GeometryInfo<dim>::face_indices())
    if (!cell->face(f)->at_boundary())
      {
        fe_iv.reinit(cell,
                     f,
                     numbers::invalid_unsigned_int,
                     cell->neighbor(f),
                     cell->neighbor_of_neighbor(f),
                     numbers::invalid_unsigned_int);
        const auto        &q_points    = fe_iv.get_quadrature_points();
        const unsigned int n_q_points  = q_points.size();
        const unsigned int n_dofs_face = fe_iv.n_current_interface_dofs();

        for (unsigned int c = 0; c < fe.n_components(); ++c)
          {
            // use a vector extractor if there
            // are sufficiently many components
            // left after the current component
            // 'c'
            if (c + dim <= fe.n_components())
              {
                const FEValuesExtractors::Vector vec_components(c);
                constexpr double                 tolerance = 1e-10;
                for (unsigned int i = 0; i < n_dofs_face; ++i)
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                      if (fe_iv[vec_components].value(true, i, q).norm() >
                          tolerance)
                        deallog << fe_iv[vec_components].value(true, i, q)
                                << "  ";
                      if (fe_iv[vec_components].value(false, i, q).norm() >
                          tolerance)
                        deallog << fe_iv[vec_components].value(false, i, q)
                                << "  ";
                      if (fe_iv[vec_components].jump_in_values(i, q).norm() >
                          tolerance)
                        deallog << fe_iv[vec_components].jump_in_values(i, q)
                                << "  ";
                      if (fe_iv[vec_components].average_of_values(i, q).norm() >
                          tolerance)
                        deallog << fe_iv[vec_components].average_of_values(i, q)
                                << "  ";
                      if (fe_iv[vec_components].jump_in_gradients(i, q).norm() >
                          tolerance)
                        deallog << fe_iv[vec_components].jump_in_gradients(i, q)
                                << "  ";
                      if (fe_iv[vec_components]
                            .average_of_gradients(i, q)
                            .norm() > tolerance)
                        deallog
                          << fe_iv[vec_components].average_of_gradients(i, q)
                          << std::endl;

                      for (unsigned int d = 0; d < dim; ++d)
                        {
                          Assert(fe_iv[vec_components].value(true, i, q)[d] ==
                                   fe_iv.shape_value(true, i, q, c + d),
                                 ExcInternalError());
                          Assert(fe_iv[vec_components].value(false, i, q)[d] ==
                                   fe_iv.shape_value(false, i, q, c + d),
                                 ExcInternalError());
                          Assert(fe_iv[vec_components].jump_in_values(i,
                                                                      q)[d] ==
                                   fe_iv.jump_in_shape_values(i, q, c + d),
                                 ExcInternalError());
                          Assert(
                            fe_iv[vec_components].average_of_values(i, q)[d] ==
                              fe_iv.average_of_shape_values(i, q, c + d),
                            ExcInternalError());
                          Assert(fe_iv[vec_components].average_of_hessians(
                                   i, q)[d] ==
                                   fe_iv.average_of_shape_hessians(i, q, c + d),
                                 ExcInternalError());
                          Assert(fe_iv[vec_components].jump_in_hessians(i,
                                                                        q)[d] ==
                                   fe_iv.jump_in_shape_hessians(i, q, c + d),
                                 ExcInternalError());

                          Assert(
                            fe_iv[vec_components].jump_in_third_derivatives(
                              i, q)[d] ==
                              fe_iv.jump_in_shape_3rd_derivatives(i, q, c + d),
                            ExcInternalError());
                        }
                    }
              }
          }
        return;
      } // at boundary
}



template <int dim>
void
test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const SphericalManifold<dim> boundary;
  tr.set_manifold(0, boundary);

  FESystem<dim> fe(FE_Q<dim>(1),
                   1,
                   FE_DGQArbitraryNodes<dim>(QGauss<1>(2)),
                   dim);
  test(tr, fe);
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}
