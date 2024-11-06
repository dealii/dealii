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


// Solve Poisson problem problem on a simplex mesh with DG and MatrixFree and
// dim components.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "./simplex_grids.h"



double
get_penalty_parameter(const unsigned int i,
                      const unsigned int j,
                      const unsigned int degree)
{
  if (degree == 1)
    {
      if (i != j)
        return 8.0;
      if (i == 0)
        return 8.0;
      if (i == 1)
        return 64.0;
    }
  else if (degree == 2)
    {
      if (i != j)
        return 32.0;
      if (i == 0)
        return 32.0;
      if (i == 1)
        return 64.0;
    }

  DEAL_II_NOT_IMPLEMENTED();

  return 0.0;
}



template <int dim, int n_components = dim>
class PoissonOperator
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<double>;
  using number     = double;

  using FECellIntegrator = FEEvaluation<dim, -1, 0, n_components, number>;
  using FEFaceIntegrator = FEFaceEvaluation<dim, -1, 0, n_components, number>;

  PoissonOperator(const MatrixFree<dim, double> &matrix_free,
                  const unsigned int             degree)
    : matrix_free(matrix_free)
    , degree(degree)
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
      [&](const auto &data, auto &dst, const auto &, const auto range) {
        FECellIntegrator phi(matrix_free, range);
        Tensor<1, dim>   value;
        for (unsigned int d = 0; d < dim; d++)
          value[d] = 0.0;
        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            phi.reinit(cell);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_value(value, q);

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
    matrix_free.template loop<VectorType, VectorType>(
      [&](const auto &data, auto &dst, const auto &src, const auto range) {
        FECellIntegrator phi(matrix_free, range);

        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            phi.reinit(cell);
            phi.gather_evaluate(src, EvaluationFlags::gradients);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            phi.integrate_scatter(EvaluationFlags::gradients, dst);
          }
      },
      [&](const auto &data, auto &dst, const auto &src, const auto range) {
        FEFaceIntegrator fe_eval(data, range, true);
        FEFaceIntegrator fe_eval_neighbor(data, range, false);

        for (unsigned int face = range.first; face < range.second; ++face)
          {
            fe_eval.reinit(face);
            fe_eval_neighbor.reinit(face);

            fe_eval.gather_evaluate(src,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);
            fe_eval_neighbor.gather_evaluate(src,
                                             EvaluationFlags::values |
                                               EvaluationFlags::gradients);
            VectorizedArray<number> sigmaF =
              get_penalty_parameter(data.get_face_active_fe_index(range, true),
                                    data.get_face_active_fe_index(range, false),
                                    degree);

            for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
              {
                Tensor<1, dim, VectorizedArray<number>> average_value =
                  (fe_eval.get_value(q) - fe_eval_neighbor.get_value(q)) * 0.5;
                Tensor<1, dim, VectorizedArray<number>> average_valgrad =
                  fe_eval.get_normal_derivative(q) +
                  fe_eval_neighbor.get_normal_derivative(q);
                average_valgrad =
                  average_value * 2. * sigmaF - average_valgrad * 0.5;
                fe_eval.submit_normal_derivative(-average_value, q);
                fe_eval_neighbor.submit_normal_derivative(-average_value, q);
                fe_eval.submit_value(average_valgrad, q);
                fe_eval_neighbor.submit_value(-average_valgrad, q);
              }
            fe_eval.integrate_scatter(EvaluationFlags::values |
                                        EvaluationFlags::gradients,
                                      dst);
            fe_eval_neighbor.integrate_scatter(EvaluationFlags::values |
                                                 EvaluationFlags::gradients,
                                               dst);
          }
      },
      [&](const auto &data, auto &dst, const auto &src, const auto range) {
        FEFaceIntegrator fe_eval(data, range, true);

        for (unsigned int face = range.first; face < range.second; ++face)
          {
            fe_eval.reinit(face);
            fe_eval.gather_evaluate(src,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);
            VectorizedArray<number> sigmaF =
              get_penalty_parameter(data.get_face_active_fe_index(range),
                                    data.get_face_active_fe_index(range),
                                    degree);

            for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
              {
                Tensor<1, dim, VectorizedArray<number>> average_value =
                  fe_eval.get_value(q);
                Tensor<1, dim, VectorizedArray<number>> average_valgrad =
                  -fe_eval.get_normal_derivative(q);
                average_valgrad += average_value * sigmaF;
                fe_eval.submit_normal_derivative(-average_value, q);
                fe_eval.submit_value(average_valgrad, q);
              }

            fe_eval.integrate_scatter(EvaluationFlags::values |
                                        EvaluationFlags::gradients,
                                      dst);
          }
      },
      dst,
      src,
      true);
  }

private:
  const MatrixFree<dim, double> &matrix_free;
  const unsigned int             degree;
};

template <int dim>
void
test(const unsigned int degree)
{
  Triangulation<dim> tria;

  unsigned int subdivisions = degree == 1 ? 16 : 8;


  GridGenerator::subdivided_hyper_cube_with_simplices(tria, subdivisions);

  FESystem<dim>      fe(FE_SimplexDGP<dim>(degree), dim);
  QGaussSimplex<dim> quad(degree + 1);
  MappingFE<dim>     mapping(FE_SimplexP<dim>(1));


  DoFHandler<dim> dof_handler(tria);

  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  constraints.close();

  const auto solve_and_postprocess =
    [&](const auto &poisson_operator,
        auto       &x,
        auto       &b) -> std::pair<unsigned int, double> {
    ReductionControl reduction_control(3000, 1e-7, 1e-2);
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
    data_out.build_patches(mapping, degree);
    data_out.write_vtu_with_pvtu_record("./",
                                        "result-" + std::to_string(dim) + "-" +
                                          std::to_string(degree) + "-" +
                                          std::to_string(version),
                                        0,
                                        MPI_COMM_WORLD);
#endif

    Vector<double> difference(tria.n_active_cells());

    deallog << "dim=" << dim << ' ';
    deallog << "degree=" << degree << ' ';

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      x,
                                      Functions::ZeroFunction<dim>(dim),
                                      difference,
                                      quad,
                                      VectorTools::L2_norm);

    const double error =
      VectorTools::compute_global_error(tria, difference, VectorTools::L2_norm);

    if (error < 0.042)
      deallog << "OK" << std::endl;
    else
      deallog << "FAIL" << std::endl;

    return {reduction_control.last_step(), reduction_control.last_value()};
  };

  const auto mf_algo = [&]() {
    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_gradients | update_values;
    additional_data.mapping_update_flags_inner_faces =
      update_gradients | update_values;
    additional_data.mapping_update_flags_boundary_faces =
      update_gradients | update_values;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim, double>::AdditionalData::none;

    MatrixFree<dim, double> matrix_free;
    matrix_free.reinit(
      mapping, dof_handler, constraints, quad, additional_data);

    PoissonOperator<dim> poisson_operator(matrix_free, degree);

    LinearAlgebra::distributed::Vector<double> x, b;
    poisson_operator.initialize_dof_vector(x);
    poisson_operator.initialize_dof_vector(b);

    poisson_operator.rhs(b);

    return solve_and_postprocess(poisson_operator, x, b);
  };

  mf_algo();
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.depth_file(1);

  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  test<2>(1);
  test<3>(1);
}
