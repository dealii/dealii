// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test FEEvaluation for assembling the Laplace matrix. Similar to
// assemble_matrix_01 but doing the whole thing in parallel with WorkStream


#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include "../tests.h"



namespace Assembly
{
  namespace Scratch
  {
    template <int dim, int fe_degree>
    struct Data
    {
      Data(const FiniteElement<dim> &fe)
        : fe_values(fe,
                    QGauss<dim>(fe.degree + 1),
                    update_values | update_gradients | update_JxW_values)
        , fe_eval(1,
                  FEEvaluation<dim, fe_degree>(fe,
                                               QGauss<1>(fe.degree + 1),
                                               fe_values.get_update_flags()))
        , cell_matrix(fe.dofs_per_cell, fe.dofs_per_cell)
        , test_matrix(cell_matrix)
      {}

      Data(const Data &data)
        : fe_values(data.fe_values.get_mapping(),
                    data.fe_values.get_fe(),
                    data.fe_values.get_quadrature(),
                    data.fe_values.get_update_flags())
        , fe_eval(data.fe_eval)
        , cell_matrix(data.cell_matrix)
        , test_matrix(data.test_matrix)
      {}

      FEValues<dim> fe_values;

      // Need to put FEEvaluation into an aligned vector because of member
      // variable alignment requirements
      AlignedVector<FEEvaluation<dim, fe_degree>> fe_eval;

      FullMatrix<double> cell_matrix, test_matrix;
    };
  } // namespace Scratch
} // namespace Assembly



// compute matrix with (\nabla v, \nabla u) + (v, 10 * u)
template <int dim, int fe_degree>
void
assemble_on_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                 Assembly::Scratch::Data<dim, fe_degree>              &data,
                 unsigned int &)
{
  const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
  const unsigned int n_q_points    = data.fe_values.get_quadrature().size();
  data.cell_matrix                 = 0;
  data.test_matrix                 = 0;
  data.fe_values.reinit(cell);

  for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          data.cell_matrix(i, j) +=
            ((data.fe_values.shape_grad(i, q_point) *
                data.fe_values.shape_grad(j, q_point) +
              10. * data.fe_values.shape_value(i, q_point) *
                data.fe_values.shape_value(j, q_point)) *
             data.fe_values.JxW(q_point));
      }

  FEEvaluation<dim, fe_degree> &fe_eval = data.fe_eval[0];
  fe_eval.reinit(cell);
  for (unsigned int i = 0; i < dofs_per_cell;
       i += VectorizedArray<double>::size())
    {
      const unsigned int n_items =
        i + VectorizedArray<double>::size() > dofs_per_cell ?
          (dofs_per_cell - i) :
          VectorizedArray<double>::size();
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        fe_eval.begin_dof_values()[j] = VectorizedArray<double>();
      for (unsigned int v = 0; v < n_items; ++v)
        fe_eval.begin_dof_values()[i + v][v] = 1.;

      fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          fe_eval.submit_value(10. * fe_eval.get_value(q), q);
          fe_eval.submit_gradient(fe_eval.get_gradient(q), q);
        }
      fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);

      for (unsigned int v = 0; v < n_items; ++v)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          data.test_matrix(fe_eval.get_internal_dof_numbering()[j],
                           fe_eval.get_internal_dof_numbering()[i + v]) =
            fe_eval.begin_dof_values()[j][v];
    }
  data.test_matrix.add(-1., data.cell_matrix);
  AssertThrow(data.test_matrix.frobenius_norm() < 1e-10, ExcInternalError());
}


void
copy_data_local_to_global(const unsigned int &)
{}


template <int dim, int fe_degree>
void
do_test(const DoFHandler<dim> &dof)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  unsigned int dummy = 0;
  WorkStream::run(dof.begin_active(),
                  dof.end(),
                  &assemble_on_cell<dim, fe_degree>,
                  &copy_data_local_to_global,
                  Assembly::Scratch::Data<dim, fe_degree>(dof.get_fe()),
                  dummy,
                  2 * MultithreadInfo::n_threads(),
                  1);
  deallog << "OK" << std::endl;
}



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_ball(tria);

  tria.refine_global(5 - dim);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  do_test<dim, fe_degree>(dof);
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
    test<2, 4>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}
