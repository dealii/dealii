// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this tests whether FEEvaluationBase::read_dof_values and
// FEEvaluationBase::read_dof_values_plain get the same data on a vector where
// constraints are distributed

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename Number   = double>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  // make function virtual to allow derived
  // classes to define a different function
  virtual void
  operator()(const MatrixFree<dim, Number> &data,
             Vector<Number> &,
             const Vector<Number>                        &src,
             const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval(data);
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, Number> fe_eval_plain(data);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(src);

        fe_eval_plain.reinit(cell);
        fe_eval_plain.read_dof_values_plain(src);

        for (unsigned int i = 0; i < fe_eval.dofs_per_cell; ++i)
          for (unsigned int j = 0; j < VectorizedArray<Number>::size(); ++j)
            {
              error += std::fabs(fe_eval.get_dof_value(i)[j] -
                                 fe_eval_plain.get_dof_value(i)[j]);
              total += std::fabs(fe_eval.get_dof_value(i)[j]);
            }
      }
  }



  void
  test_functions(const Vector<Number> &src) const
  {
    error = 0;
    total = 0;
    Vector<Number> dst_dummy;
    data.cell_loop(
      &MatrixFreeTest<dim, fe_degree, n_q_points_1d, Number>::operator(),
      this,
      dst_dummy,
      src);

    deallog << "Error read_dof_values vs read_dof_values_plain: "
            << error / total << std::endl
            << std::endl;
  };

protected:
  const MatrixFree<dim, Number> &data;
  mutable double                 error, total;
};



template <int dim, int fe_degree, typename number>
void
do_test(const DoFHandler<dim>           &dof,
        const AffineConstraints<double> &constraints)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells()
  //          << std::endl;
  // std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
  // std::cout << "Number of constraints: " << constraints.n_constraints() <<
  // std::endl;

  Vector<number> solution(dof.n_dofs());

  // create vector with random entries
  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      if (constraints.is_constrained(i))
        continue;
      const double entry = random_value<double>();
      solution(i)        = entry;
    }

  constraints.distribute(solution);
  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.store_plain_indices   = true;
    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, fe_degree, fe_degree + 1, number> mf(mf_data);
  mf.test_functions(solution);
}

template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), 1., 2., 96, true);

  // refine a few cells
  for (unsigned int i = 0; i < 11 - 3 * dim; ++i)
    {
      unsigned int counter = 0;
      for (const auto &cell : tria.active_cell_iterators())
        {
          if (counter % (7 - i) == 0)
            cell->set_refine_flag();
          ++counter;
        }
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           1,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  do_test<dim, fe_degree, double>(dof, constraints);
}


int
main()
{
  initlog();

  deallog << std::setprecision(3);
  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    test<2, 3>();
    test<2, 4>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}
