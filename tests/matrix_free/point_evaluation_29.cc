/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Test FCL integration of FEFacePointEvaluation against old version.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/non_matching/mapping_info.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;

template <int dim, typename Number>
void
do_test(unsigned int degree)
{
  using VectorType          = LinearAlgebra::distributed::Vector<Number>;
  using VectorizedArrayType = VectorizedArray<Number>;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(FE_DGQ<dim>(degree), dim + 1));

  AffineConstraints<Number> constraints;
  constraints.close();

  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
  data.mapping_update_flags                = update_values;
  data.mapping_update_flags_inner_faces    = update_values;
  data.mapping_update_flags_boundary_faces = update_values;

  MappingQ1<dim> mapping;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(
    mapping, dof_handler, constraints, QGauss<dim>(degree + 1), data);


  // setup NM::MappingInfo
  constexpr unsigned int n_lanes = VectorizedArray<Number>::size();

  std::vector<Quadrature<dim - 1>> quadrature_vector;
  quadrature_vector.reserve((matrix_free.n_inner_face_batches() +
                             matrix_free.n_boundary_face_batches()) *
                            n_lanes);

  std::vector<std::pair<typename DoFHandler<dim>::cell_iterator, unsigned int>>
    vector_face_accessors;
  vector_face_accessors.reserve((matrix_free.n_inner_face_batches() +
                                 matrix_free.n_boundary_face_batches()) *
                                n_lanes);

  // fill container for inner face batches
  unsigned int face_batch = 0;
  for (unsigned int face_batch = 0;
       face_batch < matrix_free.n_inner_face_batches() +
                      matrix_free.n_boundary_face_batches();
       ++face_batch)
    {
      for (unsigned int v = 0; v < n_lanes; ++v)
        {
          if (v < matrix_free.n_active_entries_per_face_batch(face_batch))
            vector_face_accessors.push_back(
              matrix_free.get_face_iterator(face_batch, v));
          else
            vector_face_accessors.push_back(
              matrix_free.get_face_iterator(face_batch, 0));

          quadrature_vector.push_back(QGauss<dim - 1>(degree + 1));
        }
    }

  NonMatching::MappingInfo<dim, dim, Number> nm_mapping_info(
    mapping, update_values | update_JxW_values);
  nm_mapping_info.reinit_faces(vector_face_accessors, quadrature_vector);

  VectorType dst_1;
  matrix_free.initialize_dof_vector(dst_1);
  VectorType dst_2;
  matrix_free.initialize_dof_vector(dst_2);
  VectorType src;
  matrix_free.initialize_dof_vector(src);


  class InitialSolution : public Function<dim>
  {
  public:
    InitialSolution()
      : Function<dim>(dim + 1)
    {}

    double
    value(const Point<dim> &p, const unsigned int comp) const final
    {
      if (comp < dim)
        return p[comp] * p[comp] * (comp + 1);
      else
        return p[0] * p[0] * (comp + 1);
    }
  };

  VectorTools::interpolate(*matrix_free.get_mapping_info().mapping,
                           matrix_free.get_dof_handler(),
                           InitialSolution(),
                           src);

  const auto boundary_function_1 = [&](const auto &data,
                                       auto       &dst,
                                       const auto &src,
                                       const auto  face_range) {
    FEFacePointEvaluation<dim, dim, dim, Number> phi_pnt(
      nm_mapping_info, data.get_dof_handler().get_fe(), true, 1);
    AlignedVector<Number> buffer(data.get_dof_handler().get_fe().dofs_per_cell);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        for (unsigned int v = 0;
             v < matrix_free.n_active_entries_per_face_batch(face);
             ++v)
          {
            const auto &[cell, f] =
              matrix_free.get_face_iterator(face, v, true);

            phi_pnt.reinit(face * n_lanes + v);

            cell->get_dof_values(src, buffer.begin(), buffer.end());
            phi_pnt.evaluate(buffer, EvaluationFlags::values);

            for (unsigned int q : phi_pnt.quadrature_point_indices())
              phi_pnt.submit_value(phi_pnt.get_value(q), q);

            phi_pnt.integrate(buffer, EvaluationFlags::values);

            cell->distribute_local_to_global(buffer.begin(), buffer.end(), dst);
          }
      }
  };

  const auto boundary_function_2 =
    [&](const auto &data, auto &dst, const auto &src, const auto face_range) {
      FEFaceEvaluation<dim, -1, 0, dim, Number, VectorizedArrayType> phi(
        matrix_free, true, 0, 0, 1);
      FEFacePointEvaluation<dim, dim, dim, Number> phi_pnt(
        nm_mapping_info, data.get_dof_handler().get_fe(), true, 1);

      for (unsigned int face = face_range.first; face < face_range.second;
           ++face)
        {
          phi.reinit(face);
          phi.read_dof_values(src);
          phi.project_to_face(EvaluationFlags::values);
          for (unsigned int v = 0;
               v < matrix_free.n_active_entries_per_face_batch(face);
               ++v)
            {
              phi_pnt.reinit(face * n_lanes + v);
              phi_pnt.template evaluate_in_face<VectorizedArrayType::size()>(
                &phi.get_scratch_data().begin()[0][v], EvaluationFlags::values);

              for (unsigned int q : phi_pnt.quadrature_point_indices())
                phi_pnt.submit_value(phi_pnt.get_value(q), q);

              phi_pnt.template integrate_in_face<VectorizedArrayType::size()>(
                &phi.get_scratch_data().begin()[0][v], EvaluationFlags::values);
            }
          phi.collect_from_face(EvaluationFlags::values);
          phi.distribute_local_to_global(dst);
        }
    };


  matrix_free.template loop<VectorType, VectorType>(
    {}, {}, boundary_function_1, dst_1, src, true);

  matrix_free.template loop<VectorType, VectorType>(
    {}, {}, boundary_function_2, dst_2, src, true);

  Assert(std::abs(dst_1.l2_norm() - dst_2.l2_norm()) < 1.0e-12,
         ExcInternalError());
}

int
main(int argc, char *argv[])
{
  initlog();

  do_test<2, double>(2);

  deallog << "OK" << std::endl;

  return 0;
}
