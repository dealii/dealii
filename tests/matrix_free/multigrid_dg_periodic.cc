// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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



// same as multigrid_dg_sip but using periodic boundary conditions

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename number   = double>
class LaplaceOperator : public Subscriptor
{
public:
  typedef number value_type;

  LaplaceOperator(){};

  void
  initialize(const Mapping<dim> &   mapping,
             const DoFHandler<dim> &dof_handler,
             const unsigned int     level = numbers::invalid_unsigned_int)
  {
    const QGauss<1>                                  quad(n_q_points_1d);
    typename MatrixFree<dim, number>::AdditionalData addit_data;
    addit_data.tasks_parallel_scheme =
      MatrixFree<dim, number>::AdditionalData::none;
    addit_data.tasks_block_size = 3;
    addit_data.mg_level         = level;
    addit_data.mapping_update_flags_inner_faces =
      (update_gradients | update_JxW_values);
    addit_data.mapping_update_flags_boundary_faces =
      (update_gradients | update_JxW_values);
    AffineConstraints<double> constraints;
    constraints.close();

    data.reinit(mapping, dof_handler, constraints, quad, addit_data);

    compute_inverse_diagonal();
  }

  void
  vmult(LinearAlgebra::distributed::Vector<number> &      dst,
        const LinearAlgebra::distributed::Vector<number> &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

  void
  Tvmult(LinearAlgebra::distributed::Vector<number> &      dst,
         const LinearAlgebra::distributed::Vector<number> &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

  void
  Tvmult_add(LinearAlgebra::distributed::Vector<number> &      dst,
             const LinearAlgebra::distributed::Vector<number> &src) const
  {
    vmult_add(dst, src);
  }

  void
  vmult_add(LinearAlgebra::distributed::Vector<number> &      dst,
            const LinearAlgebra::distributed::Vector<number> &src) const
  {
    Assert(src.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
    Assert(dst.partitioners_are_globally_compatible(
             *data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
    data.loop(&LaplaceOperator::local_apply,
              &LaplaceOperator::local_apply_face,
              &LaplaceOperator::local_apply_boundary,
              this,
              dst,
              src);
  }

  types::global_dof_index
  m() const
  {
    return data.get_vector_partitioner()->size();
  }

  types::global_dof_index
  n() const
  {
    return data.get_vector_partitioner()->size();
  }

  number
  el(const unsigned int row, const unsigned int col) const
  {
    AssertThrow(false,
                ExcMessage("Matrix-free does not allow for entry access"));
    return number();
  }

  void
  initialize_dof_vector(
    LinearAlgebra::distributed::Vector<number> &vector) const
  {
    if (!vector.partitioners_are_compatible(
          *data.get_dof_info(0).vector_partitioner))
      data.initialize_dof_vector(vector);
  }

  const LinearAlgebra::distributed::Vector<number> &
  get_matrix_diagonal_inverse() const
  {
    return inverse_diagonal_entries;
  }


private:
  void
  local_apply(const MatrixFree<dim, number> &                   data,
              LinearAlgebra::distributed::Vector<number> &      dst,
              const LinearAlgebra::distributed::Vector<number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(false, true, false);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate(false, true);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, number> &                   data,
    LinearAlgebra::distributed::Vector<number> &      dst,
    const LinearAlgebra::distributed::Vector<number> &src,
    const std::pair<unsigned int, unsigned int> &     face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> fe_eval(data,
                                                                       true);
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> fe_eval_neighbor(
      data, false);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        fe_eval.reinit(face);
        fe_eval_neighbor.reinit(face);

        fe_eval.read_dof_values(src);
        fe_eval.evaluate(true, true);
        fe_eval_neighbor.read_dof_values(src);
        fe_eval_neighbor.evaluate(true, true);
        VectorizedArray<number> sigmaF =
          (std::abs((fe_eval.get_normal_vector(0) *
                     fe_eval.inverse_jacobian(0))[dim - 1]) +
           std::abs((fe_eval.get_normal_vector(0) *
                     fe_eval_neighbor.inverse_jacobian(0))[dim - 1])) *
          (number)(std::max(fe_degree, 1) * (fe_degree + 1.0));

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
        fe_eval.integrate(true, true);
        fe_eval.distribute_local_to_global(dst);
        fe_eval_neighbor.integrate(true, true);
        fe_eval_neighbor.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_boundary(
    const MatrixFree<dim, number> &                   data,
    LinearAlgebra::distributed::Vector<number> &      dst,
    const LinearAlgebra::distributed::Vector<number> &src,
    const std::pair<unsigned int, unsigned int> &     face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> fe_eval(data,
                                                                       true);
    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        fe_eval.reinit(face);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(true, true);
        VectorizedArray<number> sigmaF =
          std::abs((fe_eval.get_normal_vector(0) *
                    fe_eval.inverse_jacobian(0))[dim - 1]) *
          (number)(std::max(1, fe_degree) * (fe_degree + 1.0)) * 2.;

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            VectorizedArray<number> average_value = fe_eval.get_value(q);
            VectorizedArray<number> average_valgrad =
              -fe_eval.get_normal_derivative(q);
            average_valgrad += average_value * sigmaF * 2.0;
            fe_eval.submit_normal_derivative(-average_value, q);
            fe_eval.submit_value(average_valgrad, q);
          }

        fe_eval.integrate(true, true);
        fe_eval.distribute_local_to_global(dst);
      }
  }

  void
  compute_inverse_diagonal()
  {
    data.initialize_dof_vector(inverse_diagonal_entries);
    unsigned int dummy;
    data.loop(&LaplaceOperator::local_diagonal_cell,
              &LaplaceOperator::local_diagonal_face,
              &LaplaceOperator::local_diagonal_boundary,
              this,
              inverse_diagonal_entries,
              dummy);

    for (unsigned int i = 0; i < inverse_diagonal_entries.local_size(); ++i)
      if (std::abs(inverse_diagonal_entries.local_element(i)) > 1e-10)
        inverse_diagonal_entries.local_element(i) =
          1. / inverse_diagonal_entries.local_element(i);
      else
        inverse_diagonal_entries.local_element(i) = 1.;
  }

  void
  local_diagonal_cell(
    const MatrixFree<dim, number> &             data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const unsigned int &,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);

        VectorizedArray<number> local_diagonal_vector[phi.static_dofs_per_cell];
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi.begin_dof_values()[j] = VectorizedArray<number>();
            phi.begin_dof_values()[i] = 1.;
            phi.evaluate(false, true, false);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            phi.integrate(false, true);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.static_dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = local_diagonal_vector[i];
        phi.distribute_local_to_global(dst);
      }
  }

  void
  local_diagonal_face(
    const MatrixFree<dim, number> &             data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const unsigned int &,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data, true);
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi_outer(data,
                                                                         false);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        phi.reinit(face);
        phi_outer.reinit(face);

        VectorizedArray<number> local_diagonal_vector[phi.static_dofs_per_cell];
        VectorizedArray<number> sigmaF =
          (std::abs(
             (phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]) +
           std::abs((phi.get_normal_vector(0) *
                     phi_outer.inverse_jacobian(0))[dim - 1])) *
          (number)(std::max(fe_degree, 1) * (fe_degree + 1.0));

        // Compute phi part
        for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
          phi_outer.begin_dof_values()[j] = VectorizedArray<number>();
        phi_outer.evaluate(true, true);
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi.begin_dof_values()[j] = VectorizedArray<number>();
            phi.begin_dof_values()[i] = 1.;
            phi.evaluate(true, true);

            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                VectorizedArray<number> average_value =
                  (phi.get_value(q) - phi_outer.get_value(q)) * 0.5;
                VectorizedArray<number> average_valgrad =
                  phi.get_normal_derivative(q) +
                  phi_outer.get_normal_derivative(q);
                average_valgrad =
                  average_value * 2. * sigmaF - average_valgrad * 0.5;
                phi.submit_normal_derivative(-average_value, q);
                phi.submit_value(average_valgrad, q);
              }
            phi.integrate(true, true);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = local_diagonal_vector[i];
        phi.distribute_local_to_global(dst);

        // Compute phi_outer part
        for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
          phi.begin_dof_values()[j] = VectorizedArray<number>();
        phi.evaluate(true, true);
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi_outer.begin_dof_values()[j] = VectorizedArray<number>();
            phi_outer.begin_dof_values()[i] = 1.;
            phi_outer.evaluate(true, true);

            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                VectorizedArray<number> average_value =
                  (phi.get_value(q) - phi_outer.get_value(q)) * 0.5;
                VectorizedArray<number> average_valgrad =
                  phi.get_normal_derivative(q) +
                  phi_outer.get_normal_derivative(q);
                average_valgrad =
                  average_value * 2. * sigmaF - average_valgrad * 0.5;
                phi_outer.submit_normal_derivative(-average_value, q);
                phi_outer.submit_value(-average_valgrad, q);
              }
            phi_outer.integrate(true, true);
            local_diagonal_vector[i] = phi_outer.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi_outer.begin_dof_values()[i] = local_diagonal_vector[i];
        phi_outer.distribute_local_to_global(dst);
      }
  }

  void
  local_diagonal_boundary(
    const MatrixFree<dim, number> &             data,
    LinearAlgebra::distributed::Vector<number> &dst,
    const unsigned int &,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, n_q_points_1d, 1, number> phi(data);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
      {
        phi.reinit(face);

        VectorizedArray<number> local_diagonal_vector[phi.static_dofs_per_cell];
        VectorizedArray<number> sigmaF =
          std::abs(
            (phi.get_normal_vector(0) * phi.inverse_jacobian(0))[dim - 1]) *
          (number)(std::max(1, fe_degree) * (fe_degree + 1.0)) * 2.;

        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi.begin_dof_values()[j] = VectorizedArray<number>();
            phi.begin_dof_values()[i] = 1.;
            phi.evaluate(true, true);

            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                VectorizedArray<number> average_value = phi.get_value(q);
                VectorizedArray<number> average_valgrad =
                  -phi.get_normal_derivative(q);
                average_valgrad += average_value * sigmaF * 2.0;
                phi.submit_normal_derivative(-average_value, q);
                phi.submit_value(average_valgrad, q);
              }

            phi.integrate(true, true);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = local_diagonal_vector[i];
        phi.distribute_local_to_global(dst);
      }
  }


  MatrixFree<dim, number>                    data;
  LinearAlgebra::distributed::Vector<number> inverse_diagonal_entries;
};



template <typename MATRIX, typename Number>
class MGCoarseIterative
  : public MGCoarseGridBase<LinearAlgebra::distributed::Vector<Number>>
{
public:
  MGCoarseIterative()
  {}

  void
  initialize(const MATRIX &matrix)
  {
    coarse_matrix = &matrix;
  }

  virtual void
  operator()(const unsigned int,
             LinearAlgebra::distributed::Vector<double> &      dst,
             const LinearAlgebra::distributed::Vector<double> &src) const
  {
    ReductionControl solver_control(1e4, 1e-50, 1e-7, false, false);
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver_coarse(
      solver_control);
    solver_coarse.solve(*coarse_matrix, dst, src, PreconditionIdentity());
  }

  const MATRIX *coarse_matrix;
};



template <typename LAPLACEOPERATOR>
class MGTransferMF
  : public MGTransferPrebuilt<
      LinearAlgebra::distributed::Vector<typename LAPLACEOPERATOR::value_type>>
{
public:
  MGTransferMF(const MGLevelObject<LAPLACEOPERATOR> &laplace)
    : laplace_operator(laplace){};

  /**
   * Overload copy_to_mg from MGTransferPrebuilt
   */
  template <int dim, class InVector, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &         mg_dof,
             MGLevelObject<LinearAlgebra::distributed::Vector<
               typename LAPLACEOPERATOR::value_type>> &dst,
             const InVector &                          src) const
  {
    for (unsigned int level = dst.min_level(); level <= dst.max_level();
         ++level)
      laplace_operator[level].initialize_dof_vector(dst[level]);
    MGTransferPrebuilt<LinearAlgebra::distributed::Vector<
      typename LAPLACEOPERATOR::value_type>>::copy_to_mg(mg_dof, dst, src);
  }

private:
  const MGLevelObject<LAPLACEOPERATOR> &laplace_operator;
};



template <int dim, int fe_degree, int n_q_points_1d, typename number>
void
do_test(const DoFHandler<dim> &dof)
{
  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  MappingQ<dim>                                          mapping(fe_degree + 1);
  LaplaceOperator<dim, fe_degree, n_q_points_1d, number> fine_matrix;
  fine_matrix.initialize(mapping, dof);

  LinearAlgebra::distributed::Vector<number> in, sol;
  fine_matrix.initialize_dof_vector(in);
  fine_matrix.initialize_dof_vector(sol);

  in = 1.;

  // set up multigrid in analogy to step-37
  typedef LaplaceOperator<dim, fe_degree, n_q_points_1d, number>
    LevelMatrixType;

  MGLevelObject<LevelMatrixType> mg_matrices;
  mg_matrices.resize(0, dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    mg_matrices[level].initialize(mapping, dof, level);


  MGCoarseIterative<LevelMatrixType, number> mg_coarse;
  mg_coarse.initialize(mg_matrices[0]);

  typedef PreconditionChebyshev<LevelMatrixType,
                                LinearAlgebra::distributed::Vector<number>>
    SMOOTHER;
  MGSmootherPrecondition<LevelMatrixType,
                         SMOOTHER,
                         LinearAlgebra::distributed::Vector<number>>
    mg_smoother;

  MGLevelObject<typename SMOOTHER::AdditionalData> smoother_data;
  smoother_data.resize(0, dof.get_triangulation().n_global_levels() - 1);
  for (unsigned int level = 0;
       level < dof.get_triangulation().n_global_levels();
       ++level)
    {
      smoother_data[level].smoothing_range     = 20.;
      smoother_data[level].degree              = 5;
      smoother_data[level].eig_cg_n_iterations = 15;
      auto preconditioner                      = std::make_shared<
        DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>>();
      preconditioner->reinit(mg_matrices[level].get_matrix_diagonal_inverse());
      smoother_data[level].preconditioner = std::move(preconditioner);
    }
  mg_smoother.initialize(mg_matrices, smoother_data);

  MGTransferMF<LevelMatrixType> mg_transfer(mg_matrices);
  mg_transfer.build(dof);

  mg::Matrix<LinearAlgebra::distributed::Vector<double>> mg_matrix(mg_matrices);

  Multigrid<LinearAlgebra::distributed::Vector<double>> mg(
    mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
  PreconditionMG<dim,
                 LinearAlgebra::distributed::Vector<double>,
                 MGTransferMF<LevelMatrixType>>
    preconditioner(dof, mg, mg_transfer);

  {
    ReductionControl control(30, 1e-20, 1e-10);
    SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);
    solver.solve(fine_matrix, sol, in, preconditioner);
  }
}



template <int dim, int fe_degree>
void
test()
{
  for (int i = 5; i < 9 - fe_degree; ++i)
    {
      parallel::distributed::Triangulation<dim> tria(
        MPI_COMM_WORLD,
        dealii::Triangulation<dim>::none,
        parallel::distributed::Triangulation<
          dim>::construct_multigrid_hierarchy);
      GridGenerator::hyper_cube(tria);
      // set boundary ids on boundaries to the number of the face
      for (unsigned int face = 2; face < GeometryInfo<dim>::faces_per_cell;
           ++face)
        tria.begin()->face(face)->set_all_boundary_ids(face);

      std::vector<
        GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
        periodic_faces;
      for (unsigned int d = 1; d < dim; ++d)
        GridTools::collect_periodic_faces(
          tria, 2 * d, 2 * d + 1, d, periodic_faces);
      tria.add_periodicity(periodic_faces);

      tria.refine_global(i - dim);

      FE_DGQ<dim>     fe(fe_degree);
      DoFHandler<dim> dof(tria);
      dof.distribute_dofs(fe);
      dof.distribute_mg_dofs();

      do_test<dim, fe_degree, fe_degree + 1, double>(dof);
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  mpi_initlog();

  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}
