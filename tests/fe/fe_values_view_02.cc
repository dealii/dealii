// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test the FEValues views and extractor classes. these tests use a primitive
// finite element and vector extractors

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
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

  const QGauss<dim> quadrature(2);
  FEValues<dim>     fe_values(fe,
                          quadrature,
                          update_values | update_gradients | update_hessians);
  fe_values.reinit(dof.begin_active());

  for (unsigned int c = 0; c < fe.n_components(); ++c)
    // use a vector extractor if there
    // are sufficiently many components
    // left after the current component
    // 'c'
    if (c + dim <= fe.n_components())
      {
        const FEValuesExtractors::Vector vec_components(c);

        for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
          for (const auto q : fe_values.quadrature_point_indices())
            {
              deallog << "i=" << i << ", q=" << q << std::endl;
              deallog << "   " << fe_values[vec_components].value(i, q) << ' ';
              for (unsigned int k = 0; k < dim; ++k)
                deallog << fe_values[vec_components].gradient(i, q)[k] << ' ';
              deallog << fe_values[vec_components].divergence(i, q) << ' ';
              for (unsigned int k = 0; k < dim; ++k)
                for (unsigned int l = 0; l < dim; ++l)
                  deallog
                    << fe_values[vec_components].symmetric_gradient(i, q)[k][l]
                    << ' ';
              deallog << std::endl;
              for (unsigned int k = 0; k < dim; ++k)
                for (unsigned int l = 0; l < dim; ++l)
                  for (unsigned int m = 0; m < dim; ++m)
                    deallog << fe_values[vec_components].hessian(i, q)[k][l][m]
                            << std::endl;

              for (unsigned int d = 0; d < dim; ++d)
                {
                  AssertThrow(fe_values[vec_components].value(i, q)[d] ==
                                fe_values.shape_value_component(i, q, c + d),
                              ExcInternalError());

                  AssertThrow(fe_values[vec_components].gradient(i, q)[d] ==
                                fe_values.shape_grad_component(i, q, c + d),
                              ExcInternalError());

                  AssertThrow(
                    fe_values[vec_components].symmetric_gradient(i, q) ==
                      decltype(fe_values[vec_components].symmetric_gradient(i,
                                                                            q))(
                        (fe_values[vec_components].gradient(i, q) +
                         transpose(fe_values[vec_components].gradient(i, q))) /
                        2),
                    ExcInternalError());

                  AssertThrow(fe_values[vec_components].hessian(i, q)[d] ==
                                fe_values.shape_hessian_component(i, q, c + d),
                              ExcInternalError());
                }

              AssertThrow(fe_values[vec_components].divergence(i, q) ==
                            trace(fe_values[vec_components].gradient(i, q)),
                          ExcInternalError());
            }
      }
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
                   FE_Q<dim>(2),
                   2,
                   FE_DGQArbitraryNodes<dim>(QIterated<1>(QTrapezoid<1>(), 3)),
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
