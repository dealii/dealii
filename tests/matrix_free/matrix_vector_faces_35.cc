// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// similar to matrix_vector_faces_34 (matrix-free face evaluation on a case
// with different orientations and refinements from both sides, system of
// equations), but for the advection operation with FE_DGQ.

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"

#include "matrix_vector_faces_common.h"



void
generate_grid(Triangulation<3> &triangulation, int orientation)
{
  Point<3>              vertices_1[] = {Point<3>(-1., -1., -3.),
                                        Point<3>(+1., -1., -3.),
                                        Point<3>(-1., +1., -3.),
                                        Point<3>(+1., +1., -3.),
                                        Point<3>(-1., -1., -1.),
                                        Point<3>(+1., -1., -1.),
                                        Point<3>(-1., +1., -1.),
                                        Point<3>(+1., +1., -1.),
                                        Point<3>(-1., -1., +1.),
                                        Point<3>(+1., -1., +1.),
                                        Point<3>(-1., +1., +1.),
                                        Point<3>(+1., +1., +1.)};
  std::vector<Point<3>> vertices(&vertices_1[0], &vertices_1[12]);

  std::vector<CellData<3>> cells(2, CellData<3>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<3>::vertices_per_cell] = {
    0, 1, 2, 3, 4, 5, 6, 7};

  /* cell 1 */
  int cell_vertices_1[8][GeometryInfo<3>::vertices_per_cell] = {
    {4, 5, 6, 7, 8, 9, 10, 11},
    {5, 7, 4, 6, 9, 11, 8, 10},
    {9, 8, 11, 10, 5, 4, 7, 6},
    {8, 10, 9, 11, 4, 6, 5, 7},
    {7, 6, 5, 4, 11, 10, 9, 8},
    {6, 4, 7, 5, 10, 8, 11, 9},
    {10, 11, 8, 9, 6, 7, 4, 5},
    {11, 9, 10, 8, 7, 5, 6, 4}};

  for (const unsigned int j : GeometryInfo<3>::vertex_indices())
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[orientation][j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;


  triangulation.create_triangulation(vertices, cells, SubCellData());

  const auto cell = ++(triangulation.begin());
  for (const auto face_no : cell->face_indices())
    if (!cell->face(face_no)->at_boundary())
      {
        deallog << "Orientation index within MatrixFree: "
                << (!cell->face_orientation(face_no) +
                    2 * cell->face_rotation(face_no) +
                    4 * cell->face_flip(face_no))
                << std::endl;
      }
}



template <int dim,
          int fe_degree,
          int n_q_points_1d            = fe_degree + 1,
          typename number              = double,
          typename VectorType          = Vector<number>,
          int n_components             = 1,
          typename VectorizedArrayType = VectorizedArray<number>>
class MatrixFreeAdvectionTest
{
public:
  MatrixFreeAdvectionTest(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    const bool                                          zero_within_loop = true,
    const unsigned int start_vector_component                            = 0)
    : data(data)
    , zero_within_loop(zero_within_loop)
    , start_vector_component(start_vector_component)
  {
    for (unsigned int d = 0; d < dim; ++d)
      advection[d] = 0.4 + 0.12 * d;
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    if (!zero_within_loop)
      dst = 0;
    data.loop(
      &MatrixFreeAdvectionTest::local_apply,
      &MatrixFreeAdvectionTest::local_apply_face,
      &MatrixFreeAdvectionTest::local_apply_boundary_face,
      this,
      dst,
      src,
      zero_within_loop,
      MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::values,
      MatrixFree<dim, number, VectorizedArrayType>::DataAccessOnFaces::values);
  }

private:
  void
  local_apply(const MatrixFree<dim, number, VectorizedArrayType> &data,
              VectorType                                         &dst,
              const VectorType                                   &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim,
                 fe_degree,
                 n_q_points_1d,
                 n_components,
                 number,
                 VectorizedArrayType>
      phi(data, 0, 0, start_vector_component);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(EvaluationFlags::values);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(multiply_by_advection(advection,
                                                    phi.get_value(q)),
                              q);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      phi_m(data, true, 0, 0, start_vector_component);
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      phi_p(data, false, 0, 0, start_vector_component);
    using value_type =
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_m.reinit(face);
        phi_m.read_dof_values(src);
        phi_m.evaluate(EvaluationFlags::values);
        phi_p.reinit(face);
        phi_p.read_dof_values(src);
        phi_p.evaluate(EvaluationFlags::values);

        for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
          {
            value_type u_minus = phi_m.get_value(q),
                       u_plus  = phi_p.get_value(q);
            const VectorizedArrayType normal_times_advection =
              advection * phi_m.normal_vector(q);
            const value_type flux_times_normal =
              make_vectorized_array<number, VectorizedArrayType::size()>(0.5) *
              ((u_minus + u_plus) * normal_times_advection +
               std::abs(normal_times_advection) * (u_minus - u_plus));
            phi_m.submit_value(-flux_times_normal, q);
            phi_p.submit_value(flux_times_normal, q);
          }

        phi_m.integrate(EvaluationFlags::values);
        phi_m.distribute_local_to_global(dst);
        phi_p.integrate(EvaluationFlags::values);
        phi_p.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_boundary_face(
    const MatrixFree<dim, number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceEvaluation<dim,
                     fe_degree,
                     n_q_points_1d,
                     n_components,
                     number,
                     VectorizedArrayType>
      fe_eval(data, true, 0, 0, start_vector_component);
    using value_type =
      typename FEFaceEvaluation<dim,
                                fe_degree,
                                n_q_points_1d,
                                n_components,
                                number,
                                VectorizedArrayType>::value_type;
    value_type u_plus = {};
    for (unsigned int d = 0; d < n_components; ++d)
      u_plus[d] = 1.3;

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        fe_eval.reinit(face);
        fe_eval.read_dof_values(src);
        fe_eval.evaluate(EvaluationFlags::values);

        for (unsigned int q = 0; q < fe_eval.n_q_points; ++q)
          {
            value_type                u_minus = fe_eval.get_value(q);
            const VectorizedArrayType normal_times_advection =
              advection * fe_eval.normal_vector(q);
            const value_type flux_times_normal =
              make_vectorized_array<number, VectorizedArrayType::size()>(0.5) *
              ((u_minus + u_plus) * normal_times_advection +
               std::abs(normal_times_advection) * (u_minus - u_plus));
            fe_eval.submit_value(-flux_times_normal, q);
          }

        fe_eval.integrate(EvaluationFlags::values);
        fe_eval.distribute_local_to_global(dst);
      }
  }

  const MatrixFree<dim, number, VectorizedArrayType> &data;
  const bool                                          zero_within_loop;
  const unsigned int                                  start_vector_component;
  Tensor<1, dim, VectorizedArrayType>                 advection;
};



void
run_test(const unsigned int fe_degree)
{
  FESystem<3> fe2(FE_DGQ<3>(fe_degree), 3);

  for (unsigned int orientation = 0; orientation < 8; ++orientation)
    {
      deallog << "Testing orientation case " << orientation << std::endl;
      for (unsigned int refine = 0; refine < 2; ++refine)
        {
          Triangulation<3> tria;
          generate_grid(tria, orientation);
          if (refine == 0)
            {
              tria.begin()->set_refine_flag();
              deallog << "Standard cell refined, oriented cell unrefined"
                      << std::endl;
            }
          else
            {
              (++tria.begin())->set_refine_flag();
              deallog << "Standard cell unrefined, oriented cell refined"
                      << std::endl;
            }
          tria.execute_coarsening_and_refinement();

          DoFHandler<3> dof(tria);
          dof.distribute_dofs(fe2);
          AffineConstraints<double> constraints;
          constraints.close();

          deallog << "Testing " << dof.get_fe().get_name();
          deallog << std::endl;

          MappingQ<3> mapping(dof.get_fe().degree + 1);

          Vector<double> in(dof.n_dofs()), out(dof.n_dofs());
          Vector<double> out_dist(out);

          // Set random seed for reproducibility
          Testing::srand(42);
          for (unsigned int i = 0; i < dof.n_dofs(); ++i)
            {
              if (constraints.is_constrained(i))
                continue;
              const double entry = Testing::rand() / (double)RAND_MAX;
              in(i)              = entry;
            }

          MatrixFree<3, double> mf_data;
          const QGauss<1>       quad(dof.get_fe().degree + 1);
          typename MatrixFree<3, double>::AdditionalData data;
          data.tasks_parallel_scheme =
            MatrixFree<3, double>::AdditionalData::none;
          data.tasks_block_size = 3;
          data.mapping_update_flags_inner_faces =
            (update_gradients | update_JxW_values);
          data.mapping_update_flags_boundary_faces =
            (update_gradients | update_JxW_values);

          mf_data.reinit(mapping, dof, constraints, quad, data);

          {
            MatrixFreeAdvectionTest<3, -1, 0, double, Vector<double>, 3> mf(
              mf_data);
            mf.vmult(out, in);
          }

          {
            MatrixFreeAdvection<3, -1, 0, double, Vector<double>, 3> mf(mf_data,
                                                                        true);
            mf.vmult(out_dist, in);

            out_dist -= out;
            const double diff_norm = out_dist.linfty_norm() / out.linfty_norm();
            deallog << "Norm of difference gather_evaluate: " << diff_norm
                    << std::endl;
          }
        }
    }
}



template <int dim, int fe_degree>
void
test()
{
  if (dim == 2)
    return;
  else
    run_test(fe_degree);
}
