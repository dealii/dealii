// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test checks the correct setup of FERemoteEvaluationCommunicator
// via the factory functions
// compute_remote_communicator_faces_point_to_point_interpolation()
// and compute_remote_communicator_faces_nitsche_type_mortaring().

#include <deal.II/base/mpi.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/fe_remote_evaluation.h>

#include "../tests.h"

using namespace dealii;

using Number = double;

using VectorType = LinearAlgebra::distributed::Vector<Number>;

constexpr unsigned int dim = 2;

template <int dim>
struct Coordinates : public Function<dim>
{
  Coordinates()
    : Function<dim>(dim, 0.0)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const final
  {
    if (component == 0)
      return p[component];
    return p[component];
  }
};


void
test_point_to_point(const MatrixFree<dim, Number> &matrix_free)
{
  const Triangulation<dim> &tria =
    matrix_free.get_dof_handler().get_triangulation();

  std::vector<std::pair<types::boundary_id, std::function<std::vector<bool>()>>>
    non_matching_faces_marked_vertices;

  non_matching_faces_marked_vertices.emplace_back(std::make_pair(1, [&]() {
    std::vector<bool> mask(tria.n_vertices(), true);

    for (const auto &cell : tria.active_cell_iterators())
      for (auto const &f : cell->face_indices())
        if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == 1)
          for (const auto v : cell->vertex_indices())
            mask[cell->vertex_index(v)] = false;

    return mask;
  }));

  auto remote_communicator =
    Utilities::compute_remote_communicator_faces_point_to_point_interpolation(
      matrix_free, non_matching_faces_marked_vertices);

  // 1) Allocate vectors and interpolate grid coordinates in DoFs.
  VectorType src;
  VectorType dst;
  matrix_free.initialize_dof_vector(src);
  matrix_free.initialize_dof_vector(dst);

  VectorTools::interpolate(*matrix_free.get_mapping_info().mapping,
                           matrix_free.get_dof_handler(),
                           Coordinates<dim>(),
                           src);

  // 2) Construct remote evaluator...
  FERemoteEvaluation<dim, dim, VectorizedArray<Number>> phi_remote_eval(
    remote_communicator, matrix_free.get_dof_handler());

  // ...and precompute remote values.
  phi_remote_eval.gather_evaluate(src, EvaluationFlags::values);

  // 3) Access values at quadrature points
  const auto boundary_function =
    [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
      FEFaceEvaluation<dim, -1, 0, dim, Number> phi(matrix_free, true);
      auto phi_remote = phi_remote_eval.get_data_accessor();

      for (unsigned int face = face_range.first; face < face_range.second;
           ++face)
        {
          if (matrix_free.get_boundary_id(face) == 1)
            {
              phi.reinit(face);
              phi_remote.reinit(face);

              for (unsigned int q : phi.quadrature_point_indices())
                {
                  const auto error =
                    phi.quadrature_point(q) - phi_remote.get_value(q);

                  // only consider active entries in global error
                  for (unsigned int d = 0; d < dim; ++d)
                    for (unsigned int v = 0;
                         v < matrix_free.n_active_entries_per_face_batch(face);
                         ++v)
                      AssertThrow(std::abs(error[d][v]) < 1e-13,
                                  ExcMessage("Error too large."));
                }
            }
        }
    };

  matrix_free.template loop<VectorType, VectorType>(
    {}, {}, boundary_function, dst, src, true);
}

void
test_nitsche_type_mortaring(const MatrixFree<dim, Number> &matrix_free)
{
  const Triangulation<dim> &tria =
    matrix_free.get_dof_handler().get_triangulation();

  std::vector<std::pair<types::boundary_id, std::function<std::vector<bool>()>>>
    non_matching_faces_marked_vertices;

  non_matching_faces_marked_vertices.emplace_back(std::make_pair(1, [&]() {
    std::vector<bool> mask(tria.n_vertices(), true);

    for (const auto &cell : tria.active_cell_iterators())
      for (auto const &f : cell->face_indices())
        if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == 1)
          for (const auto v : cell->vertex_indices())
            mask[cell->vertex_index(v)] = false;

    return mask;
  }));

  typename NonMatching::MappingInfo<dim, dim, Number>::AdditionalData
    additional_data;
  additional_data.use_global_weights = true;

  NonMatching::MappingInfo<dim, dim, Number> nm_mapping_info(
    *matrix_free.get_mapping_info().mapping,
    update_quadrature_points,
    additional_data);

  auto remote_communicator =
    Utilities::compute_remote_communicator_faces_nitsche_type_mortaring(
      matrix_free, non_matching_faces_marked_vertices, 4, 0, &nm_mapping_info);

  // 1) Allocate vectors and interpolate grid coordinates in DoFs.
  VectorType src;
  VectorType dst;
  matrix_free.initialize_dof_vector(src);
  matrix_free.initialize_dof_vector(dst);

  VectorTools::interpolate(*matrix_free.get_mapping_info().mapping,
                           matrix_free.get_dof_handler(),
                           Coordinates<dim>(),
                           src);

  // 2) Construct remote evaluator...
  FERemoteEvaluation<dim, dim, Number> phi_remote_eval(
    remote_communicator, matrix_free.get_dof_handler());

  // ...and precompute remote values.
  phi_remote_eval.gather_evaluate(src, EvaluationFlags::values);

  // 3) Access values at quadrature points
  const auto boundary_function =
    [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
      FEFacePointEvaluation<dim, dim, dim, Number> phi(
        nm_mapping_info, matrix_free.get_dof_handler().get_fe());
      auto phi_remote = phi_remote_eval.get_data_accessor();


      for (unsigned int face = face_range.first; face < face_range.second;
           ++face)
        {
          if (matrix_free.get_boundary_id(face) == 1)
            {
              for (unsigned int v = 0;
                   v < matrix_free.n_active_entries_per_face_batch(face);
                   ++v)
                {
                  constexpr unsigned int n_lanes =
                    VectorizedArray<Number>::size();
                  phi.reinit(face * n_lanes + v);
                  phi_remote.reinit(face * n_lanes + v);

                  for (unsigned int q : phi.quadrature_point_indices())
                    {
                      const auto             point = phi.real_point(q);
                      Tensor<1, dim, Number> temp;
                      for (unsigned int i = 0; i < dim; ++i)
                        temp[i] = point[i];

                      const auto error = temp - phi_remote.get_value(q);
                      for (unsigned int d = 0; d < dim; ++d)
                        AssertThrow(std::abs(error[d]) < 1e-13,
                                    ExcMessage("Error too large."));
                    }
                }
            }
        }
    };

  matrix_free.template loop<VectorType, VectorType>(
    {}, {}, boundary_function, dst, src, true);
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  // Construct non-matching triangulation.
  const double       length = 1.0;
  Triangulation<dim> tria_left;
  const unsigned int subdiv_left = 11;
  GridGenerator::subdivided_hyper_rectangle(tria_left,
                                            {subdiv_left, 2 * subdiv_left},
                                            {0.0, 0.0},
                                            {0.525 * length, length});

  for (const auto &face : tria_left.active_face_iterators())
    if (face->at_boundary())
      {
        face->set_boundary_id(0);
        if (face->center()[0] > 0.525 * length - 1e-6)
          face->set_boundary_id(1);
      }

  Triangulation<dim> tria_right;
  const unsigned int subdiv_right = 4;
  GridGenerator::subdivided_hyper_rectangle(tria_right,
                                            {subdiv_right, 2 * subdiv_right},
                                            {0.525 * length, 0.0},
                                            {length, length});

  for (const auto &face : tria_right.active_face_iterators())
    if (face->at_boundary())
      face->set_boundary_id(0);

  // Merge triangulations with tolerance 0 to ensure no vertices
  // are merged.
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::merge_triangulations(tria_left,
                                      tria_right,
                                      tria,
                                      /*tolerance*/ 0.,
                                      /*copy_manifold_ids*/ false,
                                      /*copy_boundary_ids*/ true);

  // Create DoFHandler
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(FE_DGQ<dim>(3), dim));

  // Constraints are not considered in this test
  AffineConstraints<Number> constraints;
  constraints.close();

  // Setup MatrixFree
  typename MatrixFree<dim, Number>::AdditionalData data;
  data.mapping_update_flags = update_values;
  // We are not considering inner faces in this test. Therefore,
  // data.mapping_update_flags_inner_faces is not set.
  data.mapping_update_flags_boundary_faces =
    update_quadrature_points | update_values;


  MatrixFree<dim, Number> matrix_free;
  matrix_free.reinit(
    MappingQ1<dim>(), dof_handler, constraints, QGauss<dim>(4), data);

  // 2) test point-to-point
  test_point_to_point(matrix_free);

  deallog << "OK" << std::endl;

  // 3) test Nitsche-type
  test_nitsche_type_mortaring(matrix_free);

  deallog << "OK" << std::endl;

  return 0;
}
