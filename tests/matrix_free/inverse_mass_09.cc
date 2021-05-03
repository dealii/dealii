// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// Tests CellwiseInverseMassMatrix::apply() with implicit definition from
// FEEvaluationBase::JxW() for vector DG elements and implicit template
// argument fe_degree=-1 in FEEvaluation, otherwise the same as
// inverse_mass_07 and inverse_mass_02

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include "../tests.h"



template <int dim, typename Number, typename VectorType = Vector<Number>>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  void
  local_mass_operator(
    const MatrixFree<dim, Number> &              data,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, -1, 0, dim, Number> fe_eval(data);
    const unsigned int                    n_q_points = fe_eval.n_q_points;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(EvaluationFlags::values);
        for (unsigned int q = 0; q < n_q_points; ++q)
          fe_eval.submit_value(fe_eval.get_value(q), q);
        fe_eval.integrate(EvaluationFlags::values);
        fe_eval.set_dof_values(dst);
      }
  }

  void
  local_inverse_mass_operator(
    const MatrixFree<dim, Number> &              data,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, -1, 0, dim, Number> fe_eval(data);
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, -1, dim, Number>
                       mass_inv(fe_eval);
    const unsigned int n_q_points = fe_eval.n_q_points;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(src);
        mass_inv.apply(fe_eval.begin_dof_values(), fe_eval.begin_dof_values());
        fe_eval.set_dof_values(dst);
      }
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    data.cell_loop(
      &MatrixFreeTest<dim, Number, VectorType>::local_mass_operator,
      this,
      dst,
      src);
  };

  void
  apply_inverse(VectorType &dst, const VectorType &src) const
  {
    data.cell_loop(
      &MatrixFreeTest<dim, Number, VectorType>::local_inverse_mass_operator,
      this,
      dst,
      src);
  };

private:
  const MatrixFree<dim, Number> &data;
};



template <int dim, typename number>
void
do_test(const DoFHandler<dim> &dof)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1> quad(dof.get_fe().degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::partition_color;
    data.tasks_block_size = 3;
    AffineConstraints<double> constraints;

    mf_data.reinit(dof, constraints, quad, data);
  }

  MatrixFreeTest<dim, number> mf(mf_data);
  Vector<number>              in(dof.n_dofs()), inverse(dof.n_dofs()),
    reference(dof.n_dofs());

  for (unsigned int i = 0; i < dof.n_dofs(); ++i)
    {
      const double entry = random_value<double>();
      in(i)              = entry;
    }

  mf.apply_inverse(inverse, in);

  SolverControl            control(1000, 1e-12);
  SolverCG<Vector<number>> solver(control);
  solver.solve(mf, reference, in, PreconditionIdentity());

  inverse -= reference;
  const double diff_norm = inverse.linfty_norm() / reference.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}



template <int dim>
void
test(const unsigned int fe_degree)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->center().norm() < 1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  const unsigned int degree = fe_degree;
  FESystem<dim>      fe(FE_DGQ<dim>(degree), dim);
  DoFHandler<dim>    dof(tria);
  dof.distribute_dofs(fe);

  do_test<dim, double>(dof);
}



int
main()
{
  initlog();

  deallog << std::setprecision(3);
  deallog.depth_file(2);

  {
    deallog.push("2d");
    test<2>(1);
    test<2>(2);
    test<2>(4);
    test<2>(9);
    deallog.pop();
    deallog.push("3d");
    test<3>(1);
    test<3>(2);
    deallog.pop();
  }
}
