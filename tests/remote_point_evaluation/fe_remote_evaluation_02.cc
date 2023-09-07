// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

// Solve Poisson problem with DG on non-matching grid using FEFaceEvaluation as
// FEEvaluationType in FERemoteEvaluation.


#include <deal.II/base/mpi.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/matrix_free/fe_remote_evaluation.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"

using namespace dealii;

template <int dim,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class PoissonOperator
{
  using This = PoissonOperator<dim, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  using FECellIntegrator =
    FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;
  using FEFaceIntegrator =
    FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;

  using NonmatchingCommunicator =
    FERemoteEvaluationCommunicator<FEFaceIntegrator, true>;
  using FERemoteFaceIntegrator = FERemoteEvaluation<NonmatchingCommunicator>;

public:
  PoissonOperator(
    const MatrixFree<dim, Number, VectorizedArrayType>       &matrix_free,
    const std::vector<std::pair<unsigned int, unsigned int>> &face_pairs)
    : matrix_free(matrix_free)
    , phi_r(phi_r_comm, matrix_free.get_dof_handler())
    , phi_r_sigma(phi_r_comm, matrix_free.get_dof_handler().get_triangulation())
    , panalty_factor(
        compute_pentaly_factor(matrix_free.get_dof_handler().get_fe().degree,
                               1.0))
  {
    // store all boundary faces in one set
    for (const auto &face_pair : face_pairs)
      faces.insert(face_pair.first);

    FEFaceIntegrator fe_face_eval(matrix_free);
    phi_r_comm.initialize_face_pairs(face_pairs, fe_face_eval);

    compute_penalty_parameters();
  }

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
        FECellIntegrator phi(data);
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
    phi_r.gather_evaluate(src,
                          EvaluationFlags::values | EvaluationFlags::gradients);

    matrix_free.loop(&This::cell_function,
                     &This::face_function,
                     &This::boundary_function,
                     this,
                     dst,
                     src,
                     true);
  }

private:
  void
  cell_function(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                VectorType                                         &dst,
                const VectorType                                   &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FECellIntegrator phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate_scatter(EvaluationFlags::gradients, dst);
      }
  }

  void
  face_function(const MatrixFree<dim, Number, VectorizedArrayType> &data,
                VectorType                                         &dst,
                const VectorType                                   &src,
                const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceIntegrator phi_m(data, true);
    FEFaceIntegrator phi_p(data, false);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_m.reinit(face);
        phi_p.reinit(face);

        phi_m.gather_evaluate(src,
                              EvaluationFlags::values |
                                EvaluationFlags::gradients);
        phi_p.gather_evaluate(src,
                              EvaluationFlags::values |
                                EvaluationFlags::gradients);

        const auto sigma_p = phi_p.read_cell_data(array_penalty_parameter);
        const auto sigma_m = phi_m.read_cell_data(array_penalty_parameter);
        const auto sigma   = std::max(sigma_m, sigma_p) * panalty_factor;
        for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
          {
            do_face_function<true>(phi_m, phi_p, sigma, q);
          }

        phi_m.integrate_scatter(EvaluationFlags::values |
                                  EvaluationFlags::gradients,
                                dst);
        phi_p.integrate_scatter(EvaluationFlags::values |
                                  EvaluationFlags::gradients,
                                dst);
      }
  }

  // While all connected cells are treated at once at all inner faces, at
  // boundary faces and non-matching faces only the inner contributions are
  // treated.
  template <bool submit_external_value, typename OtherInt>
  void
  do_face_function(FEFaceIntegrator                            &phi_m,
                   OtherInt                                    &phi_p,
                   const typename FEFaceIntegrator::value_type &sigma,
                   const unsigned int                           q) const
  {
    const auto value_m = phi_m.get_value(q);
    const auto value_p = phi_p.get_value(q);

    const auto gradient_m = phi_m.get_gradient(q);
    const auto gradient_p = phi_p.get_gradient(q);

    const auto jump_value = (value_m - value_p) * 0.5;
    const auto avg_gradient =
      phi_m.get_normal_vector(q) * (gradient_m + gradient_p) * 0.5;

    phi_m.submit_normal_derivative(-jump_value, q);
    phi_m.submit_value(jump_value * sigma * 2.0 - avg_gradient, q);

    if constexpr (submit_external_value)
      {
        phi_p.submit_normal_derivative(-jump_value, q);
        phi_p.submit_value(-jump_value * sigma * 2.0 + avg_gradient, q);
      }
  }

  void
  boundary_function(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceIntegrator phi_m(data, true);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_m.reinit(face);
        phi_m.gather_evaluate(src,
                              EvaluationFlags::values |
                                EvaluationFlags::gradients);

        const auto sigma_m = phi_m.read_cell_data(array_penalty_parameter);

        if (is_internal_face(face))
          {
            phi_r.reinit(face);
            phi_r_sigma.reinit(face);

            for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
              {
                const auto sigma =
                  std::max(sigma_m, phi_r_sigma.get_value(q)) * panalty_factor;
                do_face_function<false>(phi_m, phi_r, sigma, q);
              }
          }
        else
          {
            // to be able to reuse do_face_function() we introduce a wrapper
            // class which returns the corret calues at boundries
            struct BoundaryIntegrator
            {
              BoundaryIntegrator(const FEFaceIntegrator &phi_m)
                : phi_m(phi_m)
              {}

              typename FEFaceIntegrator::gradient_type
              get_gradient(const unsigned int q) const
              {
                return phi_m.get_gradient(q);
              }
              typename FEFaceIntegrator::value_type
              get_value(const unsigned int q) const
              {
                return -phi_m.get_value(q);
              }

            private:
              const FEFaceIntegrator &phi_m;
            };
            BoundaryIntegrator phi_bnd(phi_m);

            for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
              {
                do_face_function<false>(phi_m, phi_bnd, sigma_m, q);
              }
          }

        phi_m.integrate_scatter(EvaluationFlags::values |
                                  EvaluationFlags::gradients,
                                dst);
      }
  }


  bool
  is_internal_face(const unsigned int face) const
  {
    return faces.find(matrix_free.get_boundary_id(face)) != faces.end();
  }

  void
  compute_penalty_parameters()
  {
    // step 1) compute penalty parameter of each cell
    const unsigned int n_cells =
      matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
    array_penalty_parameter.resize(n_cells);

    const Mapping<dim>       &mapping = *matrix_free.get_mapping_info().mapping;
    const FiniteElement<dim> &fe      = matrix_free.get_dof_handler().get_fe();
    const unsigned int        degree  = fe.degree;

    QGauss<dim>   quadrature(degree + 1);
    FEValues<dim> fe_values(mapping, fe, quadrature, update_JxW_values);

    QGauss<dim - 1>   face_quadrature(degree + 1);
    FEFaceValues<dim> fe_face_values(mapping,
                                     fe,
                                     face_quadrature,
                                     update_JxW_values);

    for (unsigned int i = 0; i < n_cells; ++i)
      for (unsigned int v = 0;
           v < matrix_free.n_active_entries_per_cell_batch(i);
           ++v)
        {
          typename DoFHandler<dim>::cell_iterator cell =
            matrix_free.get_cell_iterator(i, v);
          fe_values.reinit(cell);

          Number volume = 0;
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            volume += fe_values.JxW(q);

          Number surface_area = 0;
          for (const auto f : cell->face_indices())
            {
              fe_face_values.reinit(cell, f);
              const Number factor =
                (cell->at_boundary(f) && !cell->has_periodic_neighbor(f)) ? 1. :
                                                                            0.5;
              for (unsigned int q = 0; q < face_quadrature.size(); ++q)
                surface_area += fe_face_values.JxW(q) * factor;
            }

          array_penalty_parameter[i][v] = surface_area / volume;
        }

    // 2) convert to dof vector
    FECellIntegrator fe_eval(matrix_free);

    Vector<Number> array_penalty_parameter_dof_vector(
      matrix_free.get_dof_handler().get_triangulation().n_active_cells());
    for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
      {
        fe_eval.reinit(cell);
        const auto val = fe_eval.read_cell_data(array_penalty_parameter);

        for (unsigned int v = 0;
             v < matrix_free.n_active_entries_per_cell_batch(cell);
             ++v)
          array_penalty_parameter_dof_vector
            [matrix_free.get_cell_iterator(cell, v)->active_cell_index()] =
              val[v];
      }

    // 3) cache penalty parameters for non-matching access
    // This has to be done only once since the penalty parameters never change
    phi_r_sigma.gather_evaluate(array_penalty_parameter_dof_vector,
                                EvaluationFlags::values);
  }

  static Number
  compute_pentaly_factor(const unsigned int degree, const Number factor)
  {
    return factor * (degree + 1.0) * (degree + 1.0);
  }

  const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free;

  NonmatchingCommunicator        phi_r_comm;
  mutable FERemoteFaceIntegrator phi_r;
  mutable FERemoteFaceIntegrator phi_r_sigma;

  const double                       panalty_factor;
  AlignedVector<VectorizedArrayType> array_penalty_parameter;

  std::set<unsigned int> faces;
};

template <int dim>
void
test(const unsigned int fe_degree, const unsigned int n_global_refinements = 2)
{
  using Number              = float;
  using VectorizedArrayType = VectorizedArray<Number>;
  using VectorType          = LinearAlgebra::distributed::Vector<Number>;

  const MappingQ1<dim> mapping;
  const FE_DGQ<dim>    fe_dgq(fe_degree);
  const QGauss<dim>    quad(fe_degree + 1);

  // create non-matching grid
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  Triangulation<dim>                        tria_0, tria_1;

  GridGenerator::subdivided_hyper_rectangle(
    tria_0, {7, 7}, {0.0, 0.0}, {1.0, 1.0}, true);

  GridGenerator::subdivided_hyper_rectangle(
    tria_1, {6, 3}, {1.0, 0.0}, {3.0, 1.0}, true);

  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      face->set_boundary_id(face->boundary_id() + 2 * dim);

  GridGenerator::merge_triangulations(tria_0, tria_1, tria, 0., false, true);

  AssertDimension(tria_0.n_vertices() + tria_1.n_vertices(), tria.n_vertices());

  tria.refine_global(n_global_refinements);

  // setup face pairs which are connected
  std::vector<std::pair<unsigned int, unsigned int>> face_pairs;
  face_pairs.emplace_back(1, 2 * dim);
  face_pairs.emplace_back(2 * dim, 1);

  // create DoFHandler
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe_dgq);

  // create MatrixFree
  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
  data.mapping_update_flags =
    update_quadrature_points | update_gradients | update_values;
  data.mapping_update_flags_boundary_faces = data.mapping_update_flags;
  data.mapping_update_flags_inner_faces    = data.mapping_update_flags;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  AffineConstraints<Number>                    constraints;

  matrix_free.reinit(mapping, dof_handler, constraints, quad, data);

  PoissonOperator<dim, Number, VectorizedArrayType> op(matrix_free, face_pairs);

  VectorType rhs, solution;

  op.initialize_dof_vector(rhs);
  op.initialize_dof_vector(solution);

  op.rhs(rhs);
  rhs.zero_out_ghost_values();

  try
    {
      ReductionControl reduction_control(10000, 1e-20, 1e-2);

      // note: we need to use GMRES, since the system is non-symmetrical
      SolverGMRES<VectorType> solver(reduction_control);
      solver.solve(op, solution, rhs, PreconditionIdentity());

      solution.print(deallog.get_file_stream());
    }
  catch (const SolverControl::NoConvergence &e)
    {
      std::cout << e.what() << std::endl;
    }


  // write computed vectors to Paraview
  if (true)
    {
      DataOutBase::VtkFlags flags;
      // flags.write_higher_order_cells = true;

      DataOut<dim> data_out;
      data_out.set_flags(flags);
      data_out.add_data_vector(dof_handler, solution, "solution");

      data_out.build_patches(
        mapping,
        fe_degree + 1,
        DataOut<dim>::CurvedCellRegion::curved_inner_cells);
      data_out.write_vtu_with_pvtu_record("./",
                                          "nonmatching_poisson",
                                          0,
                                          MPI_COMM_WORLD);
    }
}

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  test<2>(2);

  return 0;
}
