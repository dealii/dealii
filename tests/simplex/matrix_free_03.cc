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

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

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


const double PENALTY = 8;


template <int dim>
class SmoothSolution : public Function<dim>
{
public:
  SmoothSolution()
    : Function<dim>()
  {}
  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<double>           &values,
             const unsigned int             component = 0) const override;
};

template <int dim>
void
SmoothSolution<dim>::value_list(const std::vector<Point<dim>> &points,
                                std::vector<double>           &values,
                                const unsigned int /*component*/) const
{
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 0.0;
}

template <int dim>
class SmoothRightHandSide : public Function<dim>
{
public:
  SmoothRightHandSide()
    : Function<dim>()
  {}
  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<double>           &values,
             const unsigned int /*component*/ = 0) const override;
};

template <int dim>
void
SmoothRightHandSide<dim>::value_list(const std::vector<Point<dim>> &points,
                                     std::vector<double>           &values,
                                     const unsigned int /*component*/) const
{
  for (unsigned int i = 0; i < values.size(); ++i)
    values[i] = 1.0;
}


template <int dim>
class PoissonOperator
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<double>;
  using number     = double;

  PoissonOperator(const MatrixFree<dim, double> &matrix_free)
    : matrix_free(matrix_free)
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
      [&](const auto &data, auto &dst, const auto &, const auto cells) {
        FEEvaluation<dim, -1, 0, 1, double> phi(data);
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

  const int fe_degree = 5; /*TODO*/


  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    matrix_free.template loop<VectorType, VectorType>(
      [&](const auto &data, auto &dst, const auto &src, const auto cell_range) {
        FEEvaluation<dim, -1, 0, 1, number> phi(data);

        for (unsigned int cell = cell_range.first; cell < cell_range.second;
             ++cell)
          {
            phi.reinit(cell);
            phi.read_dof_values(src);
            phi.evaluate(EvaluationFlags::gradients);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            phi.integrate(EvaluationFlags::gradients);
            phi.set_dof_values(dst);
          }
      },
      [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
        FEFaceEvaluation<dim, -1, 0, 1, number> fe_eval(data, face_range, true);
        FEFaceEvaluation<dim, -1, 0, 1, number> fe_eval_neighbor(data,
                                                                 face_range,
                                                                 false);

        for (unsigned int face = face_range.first; face < face_range.second;
             face++)
          {
            fe_eval.reinit(face);
            fe_eval_neighbor.reinit(face);

            fe_eval.gather_evaluate(src,
                                    EvaluationFlags::values |
                                      EvaluationFlags::gradients);
            fe_eval_neighbor.gather_evaluate(src,
                                             EvaluationFlags::values |
                                               EvaluationFlags::gradients);
            VectorizedArray<number> sigmaF = PENALTY;
            //  (std::abs((fe_eval.normal_vector(0) *
            //             fe_eval.inverse_jacobian(0))[dim - 1]) +
            //   std::abs((fe_eval.normal_vector(0) *
            //             fe_eval_neighbor.inverse_jacobian(0))[dim - 1])) *
            //  (number)(std::max(fe_degree, 1) * (fe_degree + 1.0));

            for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
              {
                VectorizedArray<number> average_value =
                  (fe_eval.get_value(q) - fe_eval_neighbor.get_value(q)) * 0.5;
                VectorizedArray<number> average_valgrad =
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
      [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
        FEFaceEvaluation<dim, -1, 0, 1, number> fe_eval(data, face_range, true);
        for (unsigned int face = face_range.first; face < face_range.second;
             face++)
          {
            fe_eval.reinit(face);
            fe_eval.read_dof_values(src);
            fe_eval.evaluate(EvaluationFlags::values |
                             EvaluationFlags::gradients);
            VectorizedArray<number> sigmaF = PENALTY;
            //  std::abs((fe_eval.normal_vector(0) *
            //            fe_eval.inverse_jacobian(0))[dim - 1]) *
            //  number(std::max(fe_degree, 1) * (fe_degree + 1.0)) * 2.0;

            for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
              {
                VectorizedArray<number> average_value = fe_eval.get_value(q);
                VectorizedArray<number> average_valgrad =
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
      src);
  }

private:
  const MatrixFree<dim, double> &matrix_free;
};


struct CopyDataFace
{
  FullMatrix<double>                   cell_matrix;
  std::vector<types::global_dof_index> joint_dof_indices;
  std::array<double, 2>                values;
  std::array<unsigned int, 2>          cell_indices;
};

struct CopyData
{
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<CopyDataFace>            face_data;
  double                               value;
  unsigned int                         cell_index;
  template <class Iterator>
  void
  reinit(const Iterator &cell, unsigned int dofs_per_cell)
  {
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
};

template <int dim>
void
test(const unsigned int degree)
{
  Triangulation<dim> tria;

#if true
  unsigned int n_subdivisions = 0;

  if (dim == 2 && degree == 1)
    n_subdivisions = 16;

  if (dim == 2 && degree == 2)
    n_subdivisions = 8;

  if (dim == 3 && degree == 1)
    n_subdivisions = 8;

  if (dim == 3 && degree == 2)
    n_subdivisions = 4;

  GridGenerator::subdivided_hyper_cube_with_simplices(tria, n_subdivisions);

  FE_SimplexDGP<dim>     fe(degree);
  QGaussSimplex<dim>     quadrature(degree + 1);
  QGaussSimplex<dim - 1> face_quadrature(degree + 1);
  MappingFE<dim>         mapping(FE_SimplexP<dim>(1));
#else
  GridGenerator::subdivided_hyper_cube(tria, dim == 2 ? 16 : 8);

  FE_DGQ<dim>    fe(degree);
  QGauss<dim>    quadrature(degree + 1);
  MappingFE<dim> mapping(FE_Q<dim>(1));

#endif

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  constraints.close();

  const auto solve_and_postprocess =
    [&](const auto &poisson_operator,
        auto       &x,
        auto       &b) -> std::pair<unsigned int, double> {
    ReductionControl reduction_control(1000, 1e-7, 1e-3);
    SolverCG<std::remove_reference_t<decltype(x)>> solver(reduction_control);

    try
      {
        solver.solve(poisson_operator, x, b, PreconditionIdentity());
      }
    catch (const std::exception &e)
      {
        deallog << e.what() << std::endl;
      }

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      printf("Solved in %d iterations.\n", reduction_control.last_step());

    constraints.distribute(x);

#if 1
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    x.update_ghost_values();
    data_out.add_data_vector(dof_handler, x, "solution");
    data_out.build_patches(mapping, 2);
    data_out.write_vtu_with_pvtu_record("./", "result", 0, MPI_COMM_WORLD);
#endif

    Vector<double> difference(tria.n_active_cells());

    deallog << "dim=" << dim << ' ';
    deallog << "degree=" << degree << ' ';

    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      x,
                                      Functions::ZeroFunction<dim>(),
                                      difference,
                                      quadrature,
                                      VectorTools::L2_norm);

    deallog << VectorTools::compute_global_error(tria,
                                                 difference,
                                                 VectorTools::L2_norm)
            << std::endl;

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
      mapping, dof_handler, constraints, quadrature, additional_data);

    PoissonOperator<dim> poisson_operator(matrix_free);

    LinearAlgebra::distributed::Vector<double> x, b;
    poisson_operator.initialize_dof_vector(x);
    poisson_operator.initialize_dof_vector(b);

    poisson_operator.rhs(b);

    return solve_and_postprocess(poisson_operator, x, b);
  };

  const auto mb_algo = [&]() {
    Vector<double> solution, system_rhs;

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    SparsityPattern        sparsity_pattern;
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    SparseMatrix<double> system_matrix;
    system_matrix.reinit(sparsity_pattern);

    const double diffusion_coefficient = 1.0;

    const auto exact_solution = std::make_shared<SmoothSolution<dim>>();
    const auto rhs_function   = std::make_shared<SmoothRightHandSide<dim>>();

    const auto cell_worker =
      [&](const auto &cell, auto &scratch_data, auto &copy_data) {
        const FEValues<dim> &fe_v          = scratch_data.reinit(cell);
        const unsigned int   dofs_per_cell = fe_v.dofs_per_cell;
        copy_data.reinit(cell, dofs_per_cell);

        const auto        &q_points    = scratch_data.get_quadrature_points();
        const unsigned int n_q_points  = q_points.size();
        const std::vector<double> &JxW = scratch_data.get_JxW_values();

        std::vector<double> rhs(n_q_points);
        rhs_function->value_list(q_points, rhs);

        for (unsigned int point = 0; point < n_q_points; ++point)
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
                copy_data.cell_matrix(i, j) +=
                  diffusion_coefficient *     // nu
                  fe_v.shape_grad(i, point) * // grad v_h
                  fe_v.shape_grad(j, point) * // grad u_h
                  JxW[point];                 // dx

              copy_data.cell_rhs(i) += fe_v.shape_value(i, point) * // v_h
                                       rhs[point] *                 // f
                                       JxW[point];                  // dx
            }
      };

    const auto boundary_worker = [&](const auto         &cell,
                                     const unsigned int &face_no,
                                     auto               &scratch_data,
                                     auto               &copy_data) {
      const FEFaceValuesBase<dim> &fe_fv = scratch_data.reinit(cell, face_no);

      const auto        &q_points      = scratch_data.get_quadrature_points();
      const unsigned int n_q_points    = q_points.size();
      const unsigned int dofs_per_cell = fe_fv.dofs_per_cell;

      const std::vector<double>         &JxW = scratch_data.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals =
        scratch_data.get_normal_vectors();

      std::vector<double> g(n_q_points);
      exact_solution->value_list(q_points, g);

      const double penalty = PENALTY;

      for (unsigned int point = 0; point < n_q_points; ++point)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              copy_data.cell_matrix(i, j) +=
                (-diffusion_coefficient *        // - nu
                   fe_fv.shape_value(i, point) * // v_h
                   (fe_fv.shape_grad(j, point) * // (grad u_h .
                    normals[point])              //  n)

                 - diffusion_coefficient *         // - nu
                     (fe_fv.shape_grad(i, point) * // (grad v_h .
                      normals[point]) *            //  n)
                     fe_fv.shape_value(j, point)   // u_h

                 + diffusion_coefficient * penalty * // + nu sigma
                     fe_fv.shape_value(i, point) *   // v_h
                     fe_fv.shape_value(j, point)     // u_h

                 ) *
                JxW[point]; // dx

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            copy_data.cell_rhs(i) +=
              (-diffusion_coefficient *        // - nu
                 (fe_fv.shape_grad(i, point) * // (grad v_h .
                  normals[point]) *            //  n)
                 g[point]                      // g


               + diffusion_coefficient * penalty *        // + nu sigma
                   fe_fv.shape_value(i, point) * g[point] // v_h g

               ) *
              JxW[point]; // dx
        }
    };

    const auto face_worker = [&](const auto         &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const auto         &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 auto               &scratch_data,
                                 auto               &copy_data) {
      const FEInterfaceValues<dim> &fe_iv =
        scratch_data.reinit(cell, f, sf, ncell, nf, nsf);

      const auto        &q_points   = fe_iv.get_quadrature_points();
      const unsigned int n_q_points = q_points.size();

      copy_data.face_data.emplace_back();
      CopyDataFace      &copy_data_face = copy_data.face_data.back();
      const unsigned int n_dofs_face    = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices  = fe_iv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);

      const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

      const double penalty = PENALTY;

      for (unsigned int point = 0; point < n_q_points; ++point)
        {
          for (unsigned int i = 0; i < n_dofs_face; ++i)
            for (unsigned int j = 0; j < n_dofs_face; ++j)
              copy_data_face.cell_matrix(i, j) +=
                (-diffusion_coefficient *                 // - nu
                   fe_iv.jump_in_shape_values(i, point) * // [v_h]
                   (fe_iv.average_of_shape_gradients(j,
                                                     point) * // ({grad u_h} .
                    normals[point])                           //  n)

                 -
                 diffusion_coefficient *                         // - nu
                   (fe_iv.average_of_shape_gradients(i, point) * // (grad v_h .
                    normals[point]) *                            //  n)
                   fe_iv.jump_in_shape_values(j, point)          // [u_h]

                 + diffusion_coefficient * penalty *        // + nu sigma
                     fe_iv.jump_in_shape_values(i, point) * // [v_h]
                     fe_iv.jump_in_shape_values(j, point)   // [u_h]

                 ) *
                JxW[point]; // dx
        }
    };

    AffineConstraints<double> constraints;
    constraints.close();
    const auto copier = [&](const auto &c) {
      constraints.distribute_local_to_global(c.cell_matrix,
                                             c.cell_rhs,
                                             c.local_dof_indices,
                                             system_matrix,
                                             system_rhs);

      for (auto &cdf : c.face_data)
        {
          constraints.distribute_local_to_global(cdf.cell_matrix,
                                                 cdf.joint_dof_indices,
                                                 system_matrix);
        }
    };

    UpdateFlags cell_flags = update_values | update_gradients |
                             update_quadrature_points | update_JxW_values;
    UpdateFlags face_flags = update_values | update_gradients |
                             update_quadrature_points | update_normal_vectors |
                             update_JxW_values;

    MeshWorker::ScratchData<dim> scratch_data(
      mapping, fe, quadrature, cell_flags, face_quadrature, face_flags);
    CopyData cd;
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          cd,
                          MeshWorker::assemble_own_cells |
                            MeshWorker::assemble_boundary_faces |
                            MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);

    return solve_and_postprocess(system_matrix, solution, system_rhs);
  };

  mb_algo();
  mf_algo();
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.depth_file(1);

  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  test<2>(/*degree=*/1);
  test<2>(/*degree=*/2);
  test<3>(/*degree=*/1);
  test<3>(/*degree=*/2);
}
