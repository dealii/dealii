// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test the evaluated values and gradients at the quadrature points for
// FE_Q_iso_Q1.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q_iso_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int n_subdivisions)
{
  const FE_Q_iso_Q1<dim> fe(n_subdivisions);
  const QIterated<1>     quad(QGauss<1>(2), n_subdivisions);

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);

  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  constraints.close();

  using number = double;
  MatrixFree<dim, number>                          matrix_free;
  typename MatrixFree<dim, number>::AdditionalData data;
  data.mapping_update_flags = update_values | update_gradients;
  matrix_free.reinit(MappingQ1<dim>(), dof, constraints, quad, data);

  FEEvaluation<dim, -1, 0, 1, number> fe_eval(matrix_free);
  fe_eval.reinit(0);

  for (unsigned int i = 0; i <= dim; ++i)
    {
      for (unsigned int j = 0; j < fe_eval.dofs_per_cell; ++j)
        {
          for (unsigned int k = 0; k < fe_eval.dofs_per_cell; ++k)
            fe_eval.begin_dof_values()[k] = static_cast<number>(j == k);

          fe_eval.evaluate(EvaluationFlags::values |
                           EvaluationFlags::gradients);

          for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
            {
              const auto temp = (i == 0) ? fe_eval.get_value(q) :
                                           fe_eval.get_gradient(q)[i - 1];
              deallog << ((std::abs(temp[0]) > 1e-8) ? 1 : 0) << " ";
            }

          deallog << std::endl;
        }

      deallog << std::endl;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  {
    deallog.push("2d");
    for (unsigned int i = 1; i <= 4; ++i)
      test<2>(i);
    deallog.pop();
    deallog.push("3d");
    for (unsigned int i = 1; i <= 4; ++i)
      test<2>(i);
    deallog.pop();
  }
}
