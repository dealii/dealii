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



// Test FEEvaluation::submit_curl() for FE_NedelecNodal: the curl-curl
// operator computed with a cell_loop() must agree with the one assembled
// with FEValues. In 2d, all cells of a batch used to get the curl of the
// first cell.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1., 1.);
  tria.refine_global(2);

  const FE_NedelecNodal<dim> fe(0);
  DoFHandler<dim>            dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  constraints.close();

  const MappingQ<dim> mapping(1);

  MatrixFree<dim, double> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraints, QGauss<1>(2));

  Vector<double> src(dof_handler.n_dofs()), dst(dof_handler.n_dofs()),
    reference(dof_handler.n_dofs());
  for (unsigned int i = 0; i < src.size(); ++i)
    src(i) = std::sin(i + 1.);

  // (curl u, curl v) with FEEvaluation:
  matrix_free.template cell_loop<Vector<double>, Vector<double>>(
    [](const auto &data, auto &dst, const auto &src, const auto &range) {
      FEEvaluation<dim, 0, 2, dim, double> phi(data);
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);
          phi.gather_evaluate(src, EvaluationFlags::gradients);
          for (const unsigned int q : phi.quadrature_point_indices())
            phi.submit_curl(phi.get_curl(q), q);
          phi.integrate_scatter(EvaluationFlags::gradients, dst);
        }
    },
    dst,
    src,
    true);

  // (curl u, curl v) with FEValues:
  FEValues<dim> fe_values(mapping,
                          fe,
                          QGauss<dim>(2),
                          update_gradients | update_JxW_values);
  const FEValuesExtractors::Vector     u(0);
  std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);
      for (const unsigned int q : fe_values.quadrature_point_indices())
        {
          typename FEValuesViews::Vector<dim>::curl_type curl_u = {};
          for (const unsigned int j : fe_values.dof_indices())
            curl_u += src(dof_indices[j]) * fe_values[u].curl(j, q);
          for (const unsigned int i : fe_values.dof_indices())
            reference(dof_indices[i]) +=
              curl_u * fe_values[u].curl(i, q) * fe_values.JxW(q);
        }
    }

  dst -= reference;
  deallog << "dim " << dim << ": relative error "
          << dst.l2_norm() / reference.l2_norm() << std::endl;
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
}
