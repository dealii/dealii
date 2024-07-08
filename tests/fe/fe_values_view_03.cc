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



// test the FEValues views and extractor classes. these tests use a
// non-primitive finite element and scalar extractors

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/vector.h>

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
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end();
       ++cell)
    {
      fe_values.reinit(cell);

      for (unsigned int c = 0; c < fe.n_components(); ++c)
        {
          const FEValuesExtractors::Scalar single_component(c);

          for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
            for (const auto q : fe_values.quadrature_point_indices())
              {
                deallog << "i=" << i << ", q=" << q << std::endl;
                deallog << "   " << fe_values[single_component].value(i, q)
                        << ' ';
                for (unsigned int k = 0; k < dim; ++k)
                  deallog << fe_values[single_component].gradient(i, q)[k]
                          << ' ';
                deallog << std::endl;
                for (unsigned int k = 0; k < dim; ++k)
                  for (unsigned int l = 0; l < dim; ++l)
                    deallog << fe_values[single_component].hessian(i, q)[k][l]
                            << std::endl;

                Assert(fe_values[single_component].value(i, q) ==
                         fe_values.shape_value_component(i, q, c),
                       ExcInternalError());

                Assert(fe_values[single_component].gradient(i, q) ==
                         fe_values.shape_grad_component(i, q, c),
                       ExcInternalError());

                Assert(fe_values[single_component].hessian(i, q) ==
                         fe_values.shape_hessian_component(i, q, c),
                       ExcInternalError());
              }
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

  FESystem<dim> fe(
    FE_Q<dim>(1), 1, FE_RaviartThomas<dim>(1), 1, FE_Nedelec<dim>(0), 1);
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
