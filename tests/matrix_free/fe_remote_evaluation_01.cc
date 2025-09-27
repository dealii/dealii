// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This test checks the access of values on boundary faces using
// FERemoteEvaluation.

#include <deal.II/base/mpi.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/fe_remote_evaluation.h>

#include "../tests.h"


using Number = double;

using VectorType = LinearAlgebra::distributed::Vector<Number>;

constexpr unsigned int dim = 2;

template <int dim>
struct Coordinates : public Function<dim>
{
  Coordinates(unsigned int offset)
    : Function<dim>(dim + offset, 0.0)
    , offset(offset)
  {}

  // Function that returns the coordinates of given point.
  // In case an offset is provided, 0.0 is returned for values
  // below the offset and the coordinates are returned for
  // components >= offset. This is used to check the functionality
  // for first selected component.
  double
  value(const Point<dim> &p, const unsigned int component) const final
  {
    if (component < offset)
      return 0.0;
    return p[component - offset];
  }

  unsigned int offset;
};

FERemoteEvaluationCommunicator<dim>
construct_comm_for_face_batches(const MatrixFree<dim, Number> &matrix_free)
{
  // Setup Communication objects for all boundary faces
  FERemoteCommunicationObjectEntityBatches<dim> co;

  // Get range of boundary face indices
  const auto face_batch_range =
    std::make_pair(matrix_free.n_inner_face_batches(),
                   matrix_free.n_inner_face_batches() +
                     matrix_free.n_boundary_face_batches());

  // Fill a quadrature vector to keep track of the quadrature sizes on each
  // face
  std::vector<unsigned int> quadrature_sizes(
    matrix_free.n_boundary_face_batches());

  // Points that are searched by rpe.
  std::vector<Point<dim>> points;

  FEFaceEvaluation<dim, -1, 0, 1, Number> phi(matrix_free, true, 0, 0, 0);
  for (unsigned int bface = 0;
       bface < face_batch_range.second - face_batch_range.first;
       ++bface)
    {
      const unsigned int face = face_batch_range.first + bface;

      phi.reinit(face);

      const unsigned int n_faces =
        matrix_free.n_active_entries_per_face_batch(face);

      co.batch_id_n_entities.push_back(std::make_pair(face, n_faces));

      // Append the quadrature points to the points we need to search
      // for.
      for (unsigned int v = 0; v < n_faces; ++v)
        {
          for (unsigned int q : phi.quadrature_point_indices())
            {
              const auto point = phi.quadrature_point(q);
              Point<dim> temp;
              for (unsigned int i = 0; i < dim; ++i)
                temp[i] = point[i][v];

              points.push_back(temp);
            }
        }

      // append quadrature size
      quadrature_sizes[bface] = phi.n_q_points;
    }

  // use rpe to search for stored points
  auto rpe =
    std::make_shared<Utilities::MPI::RemotePointEvaluation<dim>>(1.0e-9);

  rpe->reinit(points,
              matrix_free.get_dof_handler().get_triangulation(),
              *matrix_free.get_mapping_info().mapping);
  Assert(rpe->all_points_found(), ExcMessage("Not all remote points found."));
  Assert(rpe->is_map_unique(), ExcMessage("The map should be unique."));

  co.rpe = rpe;

  // Renit the communicator `FERemoteEvaluationCommunicator`
  // with the communication objects.
  FERemoteEvaluationCommunicator<dim> remote_communicator;
  remote_communicator.reinit_faces({co}, face_batch_range, quadrature_sizes);

  return remote_communicator;
}

void
test_face_batches(const MatrixFree<dim, Number> &matrix_free,
                  unsigned int                   max_first_selected_comp)
{
  // 0) Get remote communicator.
  auto remote_communicator = construct_comm_for_face_batches(matrix_free);

  for (unsigned int first_selected_comp = 0;
       first_selected_comp < max_first_selected_comp;
       ++first_selected_comp)
    {
      // 1) Allocate vectors and interpolate grid coordinates in DoFs.
      VectorType src;
      VectorType dst;
      matrix_free.initialize_dof_vector(src, first_selected_comp);
      matrix_free.initialize_dof_vector(dst, first_selected_comp);

      VectorTools::interpolate(*matrix_free.get_mapping_info().mapping,
                               matrix_free.get_dof_handler(first_selected_comp),
                               Coordinates<dim>(first_selected_comp),
                               src);

      // 2) Construct remote evaluator...
      FERemoteEvaluation<dim, dim, VectorizedArray<Number>> phi_remote_eval(
        remote_communicator,
        matrix_free.get_dof_handler(first_selected_comp),
        first_selected_comp);

      // ...and precompute remote values.
      phi_remote_eval.gather_evaluate(src, EvaluationFlags::values);

      // 3) Access values at quadrature points
      const auto boundary_function = [&](const auto &data,
                                         auto       &dst,
                                         const auto &src,
                                         const auto  face_range) {
        FEFaceEvaluation<dim, -1, 0, dim, Number> phi(
          matrix_free, true, first_selected_comp, 0, first_selected_comp);
        auto phi_remote = phi_remote_eval.get_data_accessor();

        for (unsigned int face = face_range.first; face < face_range.second;
             ++face)
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
      };

      matrix_free.template loop<VectorType, VectorType>(
        {}, {}, boundary_function, dst, src, true);
    }
}


FERemoteEvaluationCommunicator<dim>
construct_comm_for_cell_face_nos(const MatrixFree<dim, Number> &matrix_free)
{
  // Setup Communication objects for all boundary faces
  FERemoteCommunicationObjectTwoLevel<dim> co;

  // Get range of boundary face indices
  const auto face_batch_range =
    std::make_pair(matrix_free.n_inner_face_batches(),
                   matrix_free.n_inner_face_batches() +
                     matrix_free.n_boundary_face_batches());

  // Fill a vector of quadrature sizes for each cell, face pair
  std::vector<std::vector<unsigned int>> quadrature_sizes;
  for (const auto &cell : matrix_free.get_dof_handler()
                            .get_triangulation()
                            .active_cell_iterators())
    quadrature_sizes.emplace_back(std::vector<unsigned int>(cell->n_faces()));

  // Points that are searched by rpe.
  std::vector<Point<dim>> points;

  FEFaceEvaluation<dim, -1, 0, 1, Number> phi(matrix_free, true, 0, 0, 0);
  for (unsigned int bface = 0;
       bface < face_batch_range.second - face_batch_range.first;
       ++bface)
    {
      const unsigned int face = face_batch_range.first + bface;
      phi.reinit(face);

      for (unsigned int v = 0;
           v < matrix_free.n_active_entries_per_face_batch(face);
           ++v)
        {
          const auto [cell, f] = matrix_free.get_face_iterator(face, v, true);
          co.cell_face_nos.push_back(std::make_pair(cell, f));

          for (unsigned int q : phi.quadrature_point_indices())
            {
              const auto point = phi.quadrature_point(q);
              Point<dim> temp;
              for (unsigned int i = 0; i < dim; ++i)
                temp[i] = point[i][v];
              points.push_back(temp);
            }

          // append correct quadrature size
          quadrature_sizes[cell->active_cell_index()][f] = phi.n_q_points;
        }
    }

  // use rpe to search for stored points
  auto rpe =
    std::make_shared<Utilities::MPI::RemotePointEvaluation<dim>>(1.0e-9);

  rpe->reinit(points,
              matrix_free.get_dof_handler().get_triangulation(),
              *matrix_free.get_mapping_info().mapping);
  Assert(rpe->all_points_found(), ExcMessage("Not all remote points found."));
  Assert(rpe->is_map_unique(), ExcMessage("The map should be unique."));

  co.rpe = rpe;

  // Renit the communicator `FERemoteEvaluationCommunicator`
  // with the communication objects.
  FERemoteEvaluationCommunicator<dim> remote_communicator;
  remote_communicator.reinit_faces(
    {co},
    matrix_free.get_dof_handler().get_triangulation().active_cell_iterators(),
    quadrature_sizes);

  return remote_communicator;
}

void
test_cell_face_nos(const MatrixFree<dim, Number> &matrix_free,
                   unsigned int                   max_first_selected_comp)
{
  // 0) Get remote communicator.
  auto remote_communicator = construct_comm_for_cell_face_nos(matrix_free);

  for (unsigned int first_selected_comp = 0;
       first_selected_comp < max_first_selected_comp;
       ++first_selected_comp)
    {
      // 1) Allocate vectors and interpolate grid coordinates in DoFs.
      VectorType src;
      VectorType dst;
      matrix_free.initialize_dof_vector(src, first_selected_comp);
      matrix_free.initialize_dof_vector(dst, first_selected_comp);

      VectorTools::interpolate(*matrix_free.get_mapping_info().mapping,
                               matrix_free.get_dof_handler(first_selected_comp),
                               Coordinates<dim>(first_selected_comp),
                               src);

      // 2) Construct remote evaluator...
      FERemoteEvaluation<dim, dim, Number> phi_remote_eval(
        remote_communicator,
        matrix_free.get_dof_handler(first_selected_comp),
        first_selected_comp);

      // ...and precompute remote values.
      phi_remote_eval.gather_evaluate(src, EvaluationFlags::values);

      // 3) Access values at quadrature points
      const auto boundary_function = [&](const auto &data,
                                         auto       &dst,
                                         const auto &src,
                                         const auto  face_range) {
        FEFaceEvaluation<dim, -1, 0, dim, Number> phi(
          matrix_free, true, first_selected_comp, 0, first_selected_comp);
        auto phi_remote = phi_remote_eval.get_data_accessor();


        for (unsigned int face = face_range.first; face < face_range.second;
             ++face)
          {
            phi.reinit(face);


            for (unsigned int v = 0;
                 v < matrix_free.n_active_entries_per_face_batch(face);
                 ++v)
              {
                const auto [cell, f] =
                  matrix_free.get_face_iterator(face, v, true);
                phi_remote.reinit(cell->active_cell_index(), f);

                for (unsigned int q : phi.quadrature_point_indices())
                  {
                    const auto             point = phi.quadrature_point(q);
                    Tensor<1, dim, Number> temp;
                    for (unsigned int i = 0; i < dim; ++i)
                      temp[i] = point[i][v];

                    const auto error = temp - phi_remote.get_value(q);
                    for (unsigned int d = 0; d < dim; ++d)
                      AssertThrow(std::abs(error[d]) < 1e-13,
                                  ExcMessage("Error too large."));
                  }
              }
          }
      };

      matrix_free.template loop<VectorType, VectorType>(
        {}, {}, boundary_function, dst, src, true);
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  // 1) setup MatrixFree for different finite elements, i.e.
  //    FESystem(FE_DGQ, dim + first_selected_component).

  const unsigned int refinements = 2;
  const unsigned int degree      = 5;

  // Construct triangulation.
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(refinements);

  // Create DoFHandlers with different number of components to
  // check different first_selected components
  std::vector<std::unique_ptr<DoFHandler<dim>>> dof_handlers;
  unsigned int                                  max_first_selected_comp = 3;
  for (unsigned int i = 0; i < max_first_selected_comp; ++i)
    {
      dof_handlers.push_back(std::make_unique<DoFHandler<dim>>(tria));
      dof_handlers.back()->distribute_dofs(
        FESystem<dim>(FE_DGQ<dim>(degree), dim + i));
    }

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


  // Construct vector of raw pointers needed to init matrix free
  std::vector<const DoFHandler<dim> *>           dof_handlers_mf;
  std::vector<const AffineConstraints<Number> *> constraints_mf;
  for (auto &dh : dof_handlers)
    {
      dof_handlers_mf.push_back(dh.get());
      constraints_mf.push_back(&constraints);
    }

  MatrixFree<dim, Number> matrix_free;
  matrix_free.reinit(MappingQ1<dim>(),
                     dof_handlers_mf,
                     constraints_mf,
                     QGauss<dim>(degree + 1),
                     data);

  // 2) test to access remote values at points in face batches
  test_face_batches(matrix_free, max_first_selected_comp);

  // 3) test to access remote values at points on faces of cells
  test_cell_face_nos(matrix_free, max_first_selected_comp);

  deallog << "OK" << std::endl;

  return 0;
}
