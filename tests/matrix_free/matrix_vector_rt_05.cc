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

// Test the correctness of the matrix-free implementation of the
// FE_RaviartThomasNodal element by evaluating values with the rudimentary
// evaluation option of initialization via cell_iterator, i.e., the slow mode
// without MatrixFree object. The test is otherwise similar to the
// matrix_vector_rt_01 case.

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"



template <int dim>
void
test(const unsigned int degree)
{
  deallog << "Testing RT element of degree " << degree << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  FE_RaviartThomasNodal<dim> fe(degree - 1);

  DoFHandler<dim> dof(tria);

  dof.distribute_dofs(fe);

  const QGauss<1>     quad(degree + 1);
  const MappingQ<dim> mapping(degree + 1);

  // create vector with random entries
  Vector<double> solution(dof.n_dofs()), initial_condition(dof.n_dofs()),
    ref(dof.n_dofs());

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    initial_condition[i] = random_value<double>();

  // Evaluate with slow path of FEEvaluation
  FEEvaluation<dim, -1, 0, dim> fe_eval(mapping,
                                        fe,
                                        quad,
                                        update_values | update_JxW_values);
  for (const auto &cell : dof.active_cell_iterators())
    {
      fe_eval.reinit(cell);
      fe_eval.read_dof_values(initial_condition);
      fe_eval.evaluate(EvaluationFlags::values);
      for (const unsigned int q : fe_eval.quadrature_point_indices())
        fe_eval.submit_value(fe_eval.get_value(q), q);
      fe_eval.integrate(EvaluationFlags::values);
      fe_eval.distribute_local_to_global(solution);
    }

  // Evaluation with FEValues

  FEValues<dim> fe_val(mapping,
                       fe,
                       Quadrature<dim>(quad),
                       update_values | update_piola | update_JxW_values);

  const unsigned int dofs_per_cell = fe_val.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  Vector<double>                       cell_rhs(fe.dofs_per_cell);

  std::vector<Tensor<1, dim>> phi_val(fe_val.n_quadrature_points);

  const FEValuesExtractors::Vector velocities(0);

  // Assemble matrix
  for (const auto &cell : dof.active_cell_iterators())
    {
      fe_val.reinit(cell);
      fe_val[velocities].get_function_values(initial_condition, phi_val);
      cell_rhs = 0;

      for (const auto q : fe_val.quadrature_point_indices())
        {
          const double JxW = fe_val.JxW(q);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_rhs(i) += (phi_val[q] * fe_val[velocities].value(i, q)) * JxW;
        }
      cell->get_dof_indices(local_dof_indices);
      AffineConstraints<double>().distribute_local_to_global(cell_rhs,
                                                             local_dof_indices,
                                                             ref);
    }

  ref -= solution;

  const double diff_norm = ref.linfty_norm() / solution.linfty_norm();
  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}


int
main()
{
  initlog();

  {
    deallog.push("2d");
    test<2>(2);
    test<2>(3);
    deallog.pop();
    deallog.push("3d");
    test<3>(2);
    deallog.pop();
  }
}
