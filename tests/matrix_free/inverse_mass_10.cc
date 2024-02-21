// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Tests CellwiseInverseMassMatrix with random dyadic coefficients, otherwise
// the same as inverse_mass_02.cc

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include "../tests.h"



template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType = Vector<Number>>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in)
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim, Number> fe_eval(data);
    const unsigned int n_q_points = fe_eval.n_q_points;

    dyadic_coefficients.resize(n_q_points);
    inverse_dyadic_coefficients.resize(n_q_points);
    create_and_invert_dyadic_coefficients();
  };

  void
  local_mass_operator(
    const MatrixFree<dim, Number>               &data,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim, Number> fe_eval(data);
    const unsigned int n_q_points = fe_eval.n_q_points;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(EvaluationFlags::values);
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const auto value_flux =
              dyadic_coefficients[q] * fe_eval.get_value(q);
            fe_eval.submit_value(value_flux, q);
          }
        fe_eval.integrate(EvaluationFlags::values);
        fe_eval.distribute_local_to_global(dst);
      }
  }

  void
  local_inverse_mass_operator(
    const MatrixFree<dim, Number>               &data,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim, Number> fe_eval(data);
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, fe_degree, dim, Number>
                                           mass_inv(fe_eval);
    const unsigned int                     n_q_points = fe_eval.n_q_points;
    AlignedVector<VectorizedArray<Number>> inverse_JxW_values(n_q_points);
    AlignedVector<Tensor<2, dim, VectorizedArray<Number>>>
      inverse_coefficients_with_JxW(n_q_points);
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        mass_inv.fill_inverse_JxW_values(inverse_JxW_values);
        for (unsigned int q = 0; q < inverse_JxW_values.size(); ++q)
          inverse_coefficients_with_JxW[q] =
            inverse_JxW_values[q] * inverse_dyadic_coefficients[q];
        fe_eval.read_dof_values(src);
        mass_inv.apply(inverse_coefficients_with_JxW,
                       fe_eval.begin_dof_values(),
                       fe_eval.begin_dof_values());
        fe_eval.distribute_local_to_global(dst);
      }
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0;
    data.cell_loop(
      &MatrixFreeTest<dim, fe_degree, Number, VectorType>::local_mass_operator,
      this,
      dst,
      src);
  };

  void
  apply_inverse(VectorType &dst, const VectorType &src) const
  {
    dst = 0;
    data.cell_loop(&MatrixFreeTest<dim, fe_degree, Number, VectorType>::
                     local_inverse_mass_operator,
                   this,
                   dst,
                   src);
  };

private:
  void
  create_and_invert_dyadic_coefficients()
  {
    for (unsigned int q = 0; q < dyadic_coefficients.size(); ++q)
      {
        for (unsigned int d1 = 0; d1 < dim; ++d1)
          {
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              {
                dyadic_coefficients[q][d1][d2] =
                  make_vectorized_array(random_value<Number>(1.0e-2, 1.0e-1));
                // make diagonal dominant because inverse is needed
                if (d1 == d2)
                  dyadic_coefficients[q][d1][d2] += 5.0;
              }
          }
        inverse_dyadic_coefficients[q] = invert(dyadic_coefficients[q]);
      }
  }

  AlignedVector<Tensor<2, dim, VectorizedArray<Number>>> dyadic_coefficients;
  AlignedVector<Tensor<2, dim, VectorizedArray<Number>>>
                                 inverse_dyadic_coefficients;
  const MatrixFree<dim, Number> &data;
};



template <int dim, int fe_degree, typename number>
void
do_test(const DoFHandler<dim> &dof)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::partition_color;
    data.tasks_block_size = 3;
    AffineConstraints<double> constraints;

    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, fe_degree, number> mf(mf_data);
  Vector<number> in(dof.n_dofs()), inverse(dof.n_dofs()),
    reference(dof.n_dofs());

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      const double entry = random_value<double>();
      in(i)              = entry;
    }

  mf.apply_inverse(inverse, in);

  SolverControl control(1000, 1e-10);
  // GMRES solver as dyadic coefficients are not necessarily symmetric
  SolverGMRES<Vector<number>> solver(control);
  solver.solve(mf, reference, in, PreconditionIdentity());

  inverse -= reference;
  const double diff_norm = inverse.linfty_norm() / reference.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}



template <int dim, int fe_degree>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  const unsigned int degree = fe_degree;
  FESystem<dim>      fe(FE_DGQ<dim>(degree), dim);
  DoFHandler<dim>    dof(tria);
  dof.distribute_dofs(fe);

  do_test<dim, fe_degree, double>(dof);
}



int
main()
{
  initlog();

  deallog << std::setprecision(3);
  deallog.depth_file(2);

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
