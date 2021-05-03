// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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



// Tests CellwiseInverseMassMatrix on vector DG elements, similar test as
// inverse_mass_02 but using different coefficients on different components

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

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



template <int dim,
          int fe_degree,
          typename Number,
          typename VectorType = Vector<Number>>
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
    FEEvaluation<dim, fe_degree, fe_degree + 1, 3, Number> fe_eval(data);
    const unsigned int n_q_points = fe_eval.n_q_points;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(EvaluationFlags::values);
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            Tensor<1, 3, VectorizedArray<double>> val = fe_eval.get_value(q);
            val[0] *= make_vectorized_array(0.8314159);
            val[1] *= make_vectorized_array(2.3);
            val[2] *= make_vectorized_array(1.98);
            fe_eval.submit_value(val, q);
          }
        fe_eval.integrate(EvaluationFlags::values);
        fe_eval.distribute_local_to_global(dst);
      }
  }

  void
  local_inverse_mass_operator(
    const MatrixFree<dim, Number> &              data,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 3, Number> fe_eval(data);
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, fe_degree, 3, Number>
                                           mass_inv(fe_eval);
    const unsigned int                     n_q_points = fe_eval.n_q_points;
    AlignedVector<VectorizedArray<Number>> inverse_coefficients(3 * n_q_points);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        fe_eval.reinit(cell);
        mass_inv.fill_inverse_JxW_values(inverse_coefficients);
        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            inverse_coefficients[q] *= make_vectorized_array(1. / 0.8314159);
            inverse_coefficients[n_q_points + q] *=
              make_vectorized_array(1. / 2.3);
            inverse_coefficients[2 * n_q_points + q] *=
              make_vectorized_array(1. / 1.98);
          }
        fe_eval.read_dof_values(src);
        mass_inv.apply(inverse_coefficients,
                       3,
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
  const MatrixFree<dim, Number> &data;
};



template <int dim, int fe_degree, typename number>
void
do_test(const DoFHandler<dim> &dof)
{
  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MappingQ<dim>           mapping(4);
  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::partition_color;
    data.tasks_block_size = 3;
    AffineConstraints<double> constraints;

    mf_data.reinit(mapping, dof, constraints, quad, data);
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

  SolverControl            control(10000, 1e-12);
  SolverCG<Vector<number>> solver(control);
  solver.solve(mf, reference, in, PreconditionIdentity());

  inverse -= reference;
  const double diff_norm = inverse.linfty_norm() / reference.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}



template <int dim, int fe_degree>
void
test()
{
  const SphericalManifold<dim> manifold;
  Triangulation<dim>           tria;
  GridGenerator::hyper_ball(tria);
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  for (; cell != endc; ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        cell->face(f)->set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);

  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  const unsigned int degree = fe_degree;
  FESystem<dim>      fe(FE_DGQ<dim>(degree), 3);
  DoFHandler<dim>    dof(tria);
  dof.distribute_dofs(fe);

  do_test<dim, fe_degree, double>(dof);
}



int
main()
{
  initlog();

  deallog << std::setprecision(3);

  // do not output CG iteration numbers to log file because these are too
  // sensitive
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
