// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2021 by the deal.II authors
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


// Test if applying a matrix vector product with a matrix containing hessians
// on faces produces the same result with FEFaceEvaluation and FEFaceValues.
// This is checked for different combinations of EvaluationFlags, FE types and
// polynomial degrees.

#include <deal.II/base/logstream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test_hessians(const unsigned int                             degree,
              const dealii::FE_Poly<dim> &                   fe,
              const dealii::EvaluationFlags::EvaluationFlags evaluation_flags)
{
  using namespace dealii;
  using VectorizedArrayType = VectorizedArray<double>;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(4 - dim);

  unsigned int counter = 0;
  for (auto cell = tria.begin_active(); cell != tria.end(); ++cell, ++counter)
    if (cell->is_locally_owned() && counter % 3 == 0)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  MappingQGeneric<dim> mapping(1);

  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

  // FEFaceEvaluation
  typename MatrixFree<dim, double, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags_inner_faces =
    update_values | update_gradients | update_hessians;

  MatrixFree<dim, double, VectorizedArrayType> matrix_free;
  matrix_free.reinit(
    mapping, dof_handler, constraints, QGauss<1>(degree + 1), additional_data);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "Working with " << fe.get_name() << " and "
            << dof_handler.n_dofs() << " dofs" << std::endl;

  LinearAlgebra::distributed::Vector<double> src, dst, dst2;
  matrix_free.initialize_dof_vector(src);
  for (auto &v : src)
    v = random_value<double>();

  matrix_free.initialize_dof_vector(dst);
  matrix_free.initialize_dof_vector(dst2);

  // Setup FEFaceValues
  QGauss<dim - 1>             quad(degree + 1);
  FEFaceValues<dim>           fe_face_values_m(mapping,
                                     fe,
                                     quad,
                                     update_values | update_gradients |
                                       update_hessians | update_JxW_values);
  FEFaceValues<dim>           fe_regface_values_p(mapping,
                                        fe,
                                        quad,
                                        update_values | update_gradients |
                                          update_hessians | update_JxW_values);
  FESubfaceValues<dim>        fe_subface_values_p(mapping,
                                           fe,
                                           quad,
                                           update_values | update_gradients |
                                             update_hessians |
                                             update_JxW_values);
  Vector<double>              solution_values_m(fe.dofs_per_cell);
  Vector<double>              solution_values_p(fe.dofs_per_cell);
  std::vector<Tensor<2, dim>> solution_hessians(quad.size());
  std::vector<Tensor<1, dim>> solution_gradients(quad.size());
  std::vector<double>         solution_values(quad.size());

  std::vector<types::global_dof_index> dof_indices_m(fe.dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_p(fe.dofs_per_cell);
  src.update_ghost_values();
  dst2 = 0;


  matrix_free.template loop<LinearAlgebra::distributed::Vector<double>,
                            LinearAlgebra::distributed::Vector<double>>(
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      (void)matrix_free;
      (void)dst;
      (void)src;
      (void)range;
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      FEFaceEvaluation<dim, -1, 0, 1, double, VectorizedArrayType> fe_eval_m(
        matrix_free, true);
      FEFaceEvaluation<dim, -1, 0, 1, double, VectorizedArrayType> fe_eval_p(
        matrix_free, false);
      for (unsigned int face = range.first; face < range.second; ++face)
        {
          // FEFaceEvaluation
          fe_eval_m.reinit(face);
          fe_eval_m.gather_evaluate(src, evaluation_flags);
          for (unsigned int q = 0; q < fe_eval_m.n_q_points; ++q)
            {
              if (evaluation_flags & EvaluationFlags::hessians)
                fe_eval_m.submit_hessian(fe_eval_m.get_hessian(q), q);
              if (evaluation_flags & EvaluationFlags::gradients)
                fe_eval_m.submit_gradient(fe_eval_m.get_gradient(q), q);
              if (evaluation_flags & EvaluationFlags::values)
                fe_eval_m.submit_value(fe_eval_m.get_value(q), q);
            }
          fe_eval_m.integrate(evaluation_flags);
          fe_eval_m.distribute_local_to_global(dst);

          fe_eval_p.reinit(face);
          fe_eval_p.gather_evaluate(src, evaluation_flags);
          for (unsigned int q = 0; q < fe_eval_p.n_q_points; ++q)
            {
              if (evaluation_flags & EvaluationFlags::hessians)
                fe_eval_p.submit_hessian(fe_eval_p.get_hessian(q), q);
              if (evaluation_flags & EvaluationFlags::gradients)
                fe_eval_p.submit_gradient(fe_eval_p.get_gradient(q), q);
              if (evaluation_flags & EvaluationFlags::values)
                fe_eval_p.submit_value(fe_eval_p.get_value(q), q);
            }
          fe_eval_p.integrate(evaluation_flags);
          fe_eval_p.distribute_local_to_global(dst);

          // FEFaceValues
          for (unsigned int v = 0;
               v < matrix_free.n_active_entries_per_face_batch(face);
               ++v)
            {
              const auto face_accessor_inside =
                matrix_free.get_face_iterator(face, v, true);
              const auto face_accessor_outside =
                matrix_free.get_face_iterator(face, v, false);

              fe_face_values_m.reinit(face_accessor_inside.first,
                                      face_accessor_inside.second);

              const FEValuesBase<dim> *fe_face_values_p = nullptr;
              if (matrix_free.get_face_info(face).subface_index <
                  GeometryInfo<dim>::max_children_per_cell)
                {
                  fe_subface_values_p.reinit(
                    face_accessor_outside.first,
                    face_accessor_outside.second,
                    matrix_free.get_face_info(face).subface_index);
                  fe_face_values_p = &fe_subface_values_p;
                }
              else
                {
                  fe_regface_values_p.reinit(face_accessor_outside.first,
                                             face_accessor_outside.second);
                  fe_face_values_p = &fe_regface_values_p;
                }

              face_accessor_inside.first->get_dof_indices(dof_indices_m);
              face_accessor_outside.first->get_dof_indices(dof_indices_p);
              constraints.get_dof_values(src,
                                         dof_indices_m.begin(),
                                         solution_values_m.begin(),
                                         solution_values_m.end());

              constraints.get_dof_values(src,
                                         dof_indices_p.begin(),
                                         solution_values_p.begin(),
                                         solution_values_p.end());

              for (unsigned int q = 0; q < quad.size(); ++q)
                {
                  double         values = 0.;
                  Tensor<1, dim> gradients;
                  Tensor<2, dim> hessians;
                  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                    {
                      if (evaluation_flags & EvaluationFlags::hessians)
                        hessians += solution_values_m(i) *
                                    fe_face_values_m.shape_hessian(i, q);
                      if (evaluation_flags & EvaluationFlags::gradients)
                        gradients += solution_values_m(i) *
                                     fe_face_values_m.shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        values += solution_values_m(i) *
                                  fe_face_values_m.shape_value(i, q);
                    }
                  solution_hessians[q]  = hessians * fe_face_values_m.JxW(q);
                  solution_gradients[q] = gradients * fe_face_values_m.JxW(q);
                  solution_values[q]    = values * fe_face_values_m.JxW(q);
                }
              for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                {
                  double sum_hessians  = 0.;
                  double sum_gradients = 0.;
                  double sum_values    = 0.;
                  for (unsigned int q = 0; q < quad.size(); ++q)
                    {
                      if (evaluation_flags & EvaluationFlags::hessians)
                        sum_hessians += double_contract<0, 0, 1, 1>(
                          solution_hessians[q],
                          fe_face_values_m.shape_hessian(i, q));
                      if (evaluation_flags & EvaluationFlags::gradients)
                        sum_gradients += solution_gradients[q] *
                                         fe_face_values_m.shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        sum_values += solution_values[q] *
                                      fe_face_values_m.shape_value(i, q);
                    }
                  solution_values_m(i) =
                    sum_hessians + sum_gradients + sum_values;
                }
              constraints.distribute_local_to_global(solution_values_m,
                                                     dof_indices_m,
                                                     dst2);

              for (unsigned int q = 0; q < quad.size(); ++q)
                {
                  double         values = 0.;
                  Tensor<1, dim> gradients;
                  Tensor<2, dim> hessians;
                  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                    {
                      if (evaluation_flags & EvaluationFlags::hessians)
                        hessians += solution_values_p(i) *
                                    fe_face_values_p->shape_hessian(i, q);
                      if (evaluation_flags & EvaluationFlags::gradients)
                        gradients += solution_values_p(i) *
                                     fe_face_values_p->shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        values += solution_values_p(i) *
                                  fe_face_values_p->shape_value(i, q);
                    }
                  solution_hessians[q]  = hessians * fe_face_values_p->JxW(q);
                  solution_gradients[q] = gradients * fe_face_values_p->JxW(q);
                  solution_values[q]    = values * fe_face_values_p->JxW(q);
                }
              for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                {
                  double sum_hessians  = 0.;
                  double sum_gradients = 0.;
                  double sum_values    = 0.;
                  for (unsigned int q = 0; q < quad.size(); ++q)
                    {
                      if (evaluation_flags & EvaluationFlags::hessians)
                        sum_hessians += double_contract<0, 0, 1, 1>(
                          solution_hessians[q],
                          fe_face_values_p->shape_hessian(i, q));
                      if (evaluation_flags & EvaluationFlags::gradients)
                        sum_gradients += solution_gradients[q] *
                                         fe_face_values_p->shape_grad(i, q);
                      if (evaluation_flags & EvaluationFlags::values)
                        sum_values += solution_values[q] *
                                      fe_face_values_p->shape_value(i, q);
                    }
                  solution_values_p(i) =
                    sum_hessians + sum_gradients + sum_values;
                }
              constraints.distribute_local_to_global(solution_values_p,
                                                     dof_indices_p,
                                                     dst2);
            }
        }
    },
    [&](
      const auto &matrix_free, auto &dst, const auto &src, const auto &range) {
      (void)matrix_free;
      (void)dst;
      (void)src;
      (void)range;
    },
    dst,
    src,
    true);

  const unsigned int end_of_print_dst =
    dof_handler.n_dofs() > 9 ? 9 : dof_handler.n_dofs();

  deallog << "dst FEE: ";
  for (unsigned int i = 0; i < end_of_print_dst; ++i)
    deallog << dst[i] << ' ';
  deallog << std::endl;

  deallog << "dst FEV: ";
  for (unsigned int i = 0; i < end_of_print_dst; ++i)
    deallog << dst2[i] << ' ';
  deallog << std::endl;

  // compare solutions of matrix vector product
  {
    dst2 -= dst;

    double error = 0.;
    if (dst.l2_norm() > 0)
      error = dst2.l2_norm() / dst.l2_norm();
    else
      error = dst2.l2_norm();

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      deallog << "FEValues verification: " << error << std::endl << std::endl;
  }
}

int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  {
    dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::hessians;
    deallog << "test_hessians_only" << std::endl;
    for (unsigned int i = 1; i < 3; ++i)
      {
        test_hessians<1>(i, dealii::FE_Q<1>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_Q<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_Q<3>(i), evaluation_flags);
        test_hessians<1>(i, dealii::FE_DGQ<1>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_DGQ<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_DGQ<3>(i), evaluation_flags);
      }
  }

  {
    dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::values | EvaluationFlags::hessians;
    deallog << "test_hessians_with_values" << std::endl;
    for (unsigned int i = 1; i < 3; ++i)
      {
        test_hessians<1>(i, dealii::FE_Q<1>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_Q<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_Q<3>(i), evaluation_flags);
        test_hessians<1>(i, dealii::FE_DGQ<1>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_DGQ<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_DGQ<3>(i), evaluation_flags);
      }
  }

  {
    dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::gradients | EvaluationFlags::hessians;
    deallog << "test_hessians_with_gradients" << std::endl;
    for (unsigned int i = 1; i < 3; ++i)
      {
        test_hessians<1>(i, dealii::FE_Q<1>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_Q<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_Q<3>(i), evaluation_flags);
        test_hessians<1>(i, dealii::FE_DGQ<1>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_DGQ<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_DGQ<3>(i), evaluation_flags);
      }
  }

  {
    dealii::EvaluationFlags::EvaluationFlags evaluation_flags =
      EvaluationFlags::values | EvaluationFlags::gradients |
      EvaluationFlags::hessians;
    deallog << "test_hessians_with_gradients_and_values" << std::endl;
    // run the last test also for cubic polynomials
    for (unsigned int i = 1; i < 4; ++i)
      {
        test_hessians<1>(i, dealii::FE_Q<1>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_Q<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_Q<3>(i), evaluation_flags);
        test_hessians<1>(i, dealii::FE_DGQ<1>(i), evaluation_flags);
        test_hessians<2>(i, dealii::FE_DGQ<2>(i), evaluation_flags);
        test_hessians<3>(i, dealii::FE_DGQ<3>(i), evaluation_flags);
      }
  }
}
