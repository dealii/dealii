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



// Check that FEEvaluation::distribute_local_to_global and related functions
// work for FEEvaluation initialized with cell_iterator.



#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);

  FE_DGQ<dim>     fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  MappingQ1<dim> mapping;

  LinearAlgebra::distributed::Vector<double> vec;
  vec.reinit(dof_handler.n_dofs());
  FEEvaluation<dim, -1, 0, 1, double, VectorizedArray<double, 1>> evaluator(
    mapping, fe, QGauss<1>(2), update_values | update_JxW_values);

  // test distribute_local_to_global()
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      evaluator.reinit(cell);
      unsigned int count = 0;
      for (const unsigned int q : evaluator.quadrature_point_indices())
        evaluator.submit_value(1.0 + 0.1 * (count++), q);
      evaluator.integrate(EvaluationFlags::values);
      evaluator.distribute_local_to_global(vec);
    }
  deallog << "Result distribute_to_global:" << std::endl;
  vec.print(deallog.get_file_stream());
  deallog << std::endl;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      evaluator.reinit(cell);
      unsigned int count = 0;
      for (const unsigned int q : evaluator.quadrature_point_indices())
        evaluator.submit_value(1.0 + 0.1 * (count++), q);
      evaluator.integrate(EvaluationFlags::values);
      evaluator.distribute_local_to_global(vec);
    }
  deallog << "Result distribute_to_global round 2:" << std::endl;
  vec.print(deallog.get_file_stream());
  deallog << std::endl;

  // test set_dof_values()
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      evaluator.reinit(cell);
      unsigned int count = 0;
      for (const unsigned int q : evaluator.quadrature_point_indices())
        evaluator.submit_value(1.0 + 0.1 * (count++), q);
      evaluator.integrate(EvaluationFlags::values);
      evaluator.set_dof_values(vec);
    }
  deallog << "Result set_dof_values:" << std::endl;
  vec.print(deallog.get_file_stream());
  deallog << std::endl;

  // test read_dof_values()
  deallog << "Result read_dof_values:" << std::endl;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      evaluator.reinit(cell);
      evaluator.read_dof_values(vec);
      evaluator.evaluate(EvaluationFlags::values);
      for (const unsigned int q : evaluator.quadrature_point_indices())
        deallog << evaluator.get_value(q) << "  ";
      deallog << std::endl;
    }
  deallog << std::endl;

  // test read_dof_values_plain()
  deallog << "Result read_dof_values_plain:" << std::endl;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      evaluator.reinit(cell);
      evaluator.read_dof_values_plain(vec);
      evaluator.evaluate(EvaluationFlags::values);
      for (const unsigned int q : evaluator.quadrature_point_indices())
        deallog << evaluator.get_value(q) << "  ";
      deallog << std::endl;
    }
  deallog << std::endl;
}


int
main()
{
  initlog();
  test<2>();
}
