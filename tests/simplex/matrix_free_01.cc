// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Solve Poisson problem and Helmholtz problem on a simplex mesh with
// continuous elements and compare results between matrix-free and matrix-based
// implementations.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "./simplex_grids.h"



template <int dim>
class PoissonOperator
{
public:
  using VectorType = Vector<double>;

  PoissonOperator(const MatrixFree<dim, double> &matrix_free,
                  const bool                     do_helmholtz)
    : matrix_free(matrix_free)
    , do_helmholtz(do_helmholtz)
  {}

  void
  initialize_dof_vector(VectorType &vec)
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  rhs(VectorType &vec) const
  {
    const int dummy = 0;

    matrix_free.template cell_loop<VectorType, int>(
      [&](const auto &, auto &dst, const auto &, const auto cells) {
        FEEvaluation<dim, -1, 0, 1, double> phi(matrix_free);
        for (unsigned int cell = cells.first; cell < cells.second; ++cell)
          {
            phi.reinit(cell);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_value(1.0, q);

            phi.integrate_scatter(EvaluationFlags::values, dst);
          }
      },
      vec,
      dummy,
      true);
  }


  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    matrix_free.template cell_loop<VectorType, VectorType>(
      [&](const auto &, auto &dst, const auto &src, const auto cells) {
        FEEvaluation<dim, -1, 0, 1, double> phi(matrix_free);
        EvaluationFlags::EvaluationFlags    fe_eval_flags =
          EvaluationFlags::gradients;
        if (do_helmholtz)
          fe_eval_flags |= EvaluationFlags::values;
        for (unsigned int cell = cells.first; cell < cells.second; ++cell)
          {
            phi.reinit(cell);
            phi.gather_evaluate(src, fe_eval_flags);

            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                if (do_helmholtz)
                  phi.submit_value(phi.get_value(q), q);

                phi.submit_gradient(phi.get_gradient(q), q);
              }

            phi.integrate_scatter(fe_eval_flags, dst);
          }
      },
      dst,
      src,
      true);
  }

private:
  const MatrixFree<dim, double> &matrix_free;
  const bool                     do_helmholtz;
};

template <int dim>
void
test(const unsigned int v, const unsigned int degree, const bool do_helmholtz)
{
  Triangulation<dim> tria;

  std::shared_ptr<FiniteElement<dim>> fe;
  std::shared_ptr<Quadrature<dim>>    quad;
  std::shared_ptr<FiniteElement<dim>> fe_mapping;

  if (v == 0)
    {
      GridGenerator::subdivided_hyper_cube_with_simplices(tria,
                                                          dim == 2 ? 16 : 8);
      fe         = std::make_shared<FE_SimplexP<dim>>(degree);
      quad       = std::make_shared<QGaussSimplex<dim>>(degree + 1);
      fe_mapping = std::make_shared<FE_SimplexP<dim>>(1);
    }
  else if (v == 1)
    {
      GridGenerator::subdivided_hyper_cube_with_wedges(tria, dim == 2 ? 16 : 8);
      fe         = std::make_shared<FE_WedgeP<dim>>(degree);
      quad       = std::make_shared<QGaussWedge<dim>>(degree + 1);
      fe_mapping = std::make_shared<FE_WedgeP<dim>>(1);
    }
  else if (v == 2)
    {
      GridGenerator::subdivided_hyper_cube_with_pyramids(tria,
                                                         dim == 2 ? 16 : 8);
      fe         = std::make_shared<FE_PyramidP<dim>>(degree);
      quad       = std::make_shared<QGaussPyramid<dim>>(degree + 1);
      fe_mapping = std::make_shared<FE_PyramidP<dim>>(1);
    }
  else
    DEAL_II_NOT_IMPLEMENTED();

  MappingFE<dim> mapping(*fe_mapping);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(*fe);

  AffineConstraints<double> constraints;
#if false
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
#else
  DoFTools::make_zero_boundary_constraints(dof_handler, 0, constraints);
#endif
  constraints.close();

  const auto solve_and_postprocess =
    [&](const auto &poisson_operator,
        auto       &x,
        auto       &b) -> std::pair<unsigned int, double> {
    ReductionControl                               reduction_control;
    SolverCG<std::remove_reference_t<decltype(x)>> solver(reduction_control);
    solver.solve(poisson_operator, x, b, PreconditionIdentity());

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      printf("Solved in %d iterations.\n", reduction_control.last_step());

    constraints.distribute(x);

#if 0
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  x.update_ghost_values();
  data_out.add_data_vector(dof_handler, x, "solution");
  data_out.build_patches(mapping, 2);
  data_out.write_vtu_with_pvtu_record("./", "result", 0, MPI_COMM_WORLD);
#endif

    return {reduction_control.last_step(), reduction_control.last_value()};
  };

  const auto mf_algo = [&]() {
    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_gradients | update_values;

    MatrixFree<dim, double> matrix_free;
    matrix_free.reinit(
      mapping, dof_handler, constraints, *quad, additional_data);

    PoissonOperator<dim> poisson_operator(matrix_free, do_helmholtz);

    Vector<double> x, b;
    poisson_operator.initialize_dof_vector(x);
    poisson_operator.initialize_dof_vector(b);

    poisson_operator.rhs(b);

    return solve_and_postprocess(poisson_operator, x, b);
  };

  const auto mb_algo = [&]() {
    Vector<double> x, b;

    x.reinit(dof_handler.n_dofs());
    b.reinit(dof_handler.n_dofs());

    SparseMatrix<double> A;

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);

    SparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from(dsp);
    A.reinit(sparsity_pattern);

    const auto flags = update_values | update_gradients | update_JxW_values;

    FEValues<dim> fe_values(mapping, *fe, *quad, flags);

    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        fe_values.reinit(cell);

        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_rhs.reinit(dofs_per_cell);

        for (const auto q : fe_values.quadrature_point_indices())
          {
            for (const auto i : fe_values.dof_indices())
              for (const auto j : fe_values.dof_indices())
                cell_matrix(i, j) += (fe_values.shape_grad(i, q) *        //
                                        fe_values.shape_grad(j, q) +      //
                                      static_cast<double>(do_helmholtz) * //
                                        fe_values.shape_value(i, q) *     //
                                        fe_values.shape_value(j, q)) *    //
                                     fe_values.JxW(q);                    //

            for (const unsigned int i : fe_values.dof_indices())
              cell_rhs(i) += (fe_values.shape_value(i, q) * //
                              1. *                          //
                              fe_values.JxW(q));            //
          }

        local_dof_indices.resize(cell->get_fe().dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, A, b);
      }

    return solve_and_postprocess(A, x, b);
  };

  const auto compare = [&](const auto result_mf, const auto result_mb) {
    AssertDimension(result_mf.first, result_mb.first);
    Assert(std::abs(result_mf.second - result_mb.second) < 1e-8,
           ExcNotImplemented());

    deallog << "dim=" << dim << ' ';
    deallog << "degree=" << degree << ' ';
    deallog << "Type=";

    if (do_helmholtz)
      deallog << "Helmholtz"
              << " : ";
    else
      deallog << "Possion  "
              << " : ";

    deallog << "Convergence step " << result_mf.first << " value "
            << result_mf.second << '.' << std::endl;
  };

  compare(mf_algo(), mb_algo());
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.depth_file(2);

  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  for (unsigned int i = 0; i <= 2; ++i)
    {
      if (i == 0)
        deallog.push("SIMPLEX");
      else if (i == 1)
        deallog.push("WEDGE  ");
      else if (i == 2)
        deallog.push("PYRAMID");
      else
        DEAL_II_NOT_IMPLEMENTED();

      if (i == 0) // 2D makes only sense for simplex
        {
          test<2>(i, /*degree=*/1, /*do_helmholtz*/ false);
          test<2>(i, /*degree=*/1, /*do_helmholtz*/ true);
          test<2>(i, /*degree=*/2, /*do_helmholtz*/ false);
          test<2>(i, /*degree=*/2, /*do_helmholtz*/ true);
        }

      test<3>(i, /*degree=*/1, /*do_helmholtz*/ false);
      test<3>(i, /*degree=*/1, /*do_helmholtz*/ true);

      if (i !=
          2) // for pyramids no quadratic elements have been implemented yet
        {
          test<3>(i, /*degree=*/2, /*do_helmholtz*/ false);
          test<3>(i, /*degree=*/2, /*do_helmholtz*/ true);
        }

      deallog.pop();
    }
}
