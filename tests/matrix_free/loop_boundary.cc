// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this function tests the correctness of the matrix-free loop with continuous
// elements when both a cell and a boundary function are given. In an initial
// implementation, we used to miss to exchange some ghost entries upon the
// compress() functionality within the matrix-free loop

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"

template <int dim, int fe_degree>
void
do_test(const unsigned int n_refine, const bool overlap_communication)
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, 0.1, 0.96);
  tria.refine_global(n_refine);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  constraints.close();

  MatrixFree<dim>                          matrix_free;
  typename MatrixFree<dim>::AdditionalData add_data;
  add_data.mapping_update_flags = update_values | update_gradients |
                                  update_JxW_values | update_quadrature_points;
  add_data.mapping_update_flags_boundary_faces =
    update_values | update_JxW_values | update_quadrature_points;
  add_data.overlap_communication_computation = overlap_communication;
  matrix_free.reinit(MappingQ1<dim>{},
                     dof_handler,
                     constraints,
                     QGauss<1>(fe.degree + 1),
                     add_data);

  LinearAlgebra::distributed::Vector<double> in, ref, test;
  matrix_free.initialize_dof_vector(in);
  for (unsigned int i = 0; i < in.get_partitioner()->locally_owned_size(); ++i)
    in.local_element(i) = in.get_partitioner()->local_to_global(i);
  matrix_free.initialize_dof_vector(ref);
  matrix_free.initialize_dof_vector(test);

  std::function<void(const MatrixFree<dim> &,
                     LinearAlgebra::distributed::Vector<double> &,
                     const LinearAlgebra::distributed::Vector<double> &,
                     const std::pair<unsigned int, unsigned int> &)>
    cell_func = [](const MatrixFree<dim>                            &data,
                   LinearAlgebra::distributed::Vector<double>       &out,
                   const LinearAlgebra::distributed::Vector<double> &in,
                   const std::pair<unsigned int, unsigned int> &range) -> void {
    FEEvaluation<dim, fe_degree> eval(data);

    for (unsigned int cell = range.first; cell < range.second; ++cell)
      {
        eval.reinit(cell);
        eval.gather_evaluate(in, EvaluationFlags::gradients);
        for (unsigned int q = 0; q < eval.n_q_points; ++q)
          {
            eval.submit_gradient(eval.get_gradient(q), q);
            eval.submit_value(eval.quadrature_point(q).square(), q);
          }
        eval.integrate_scatter(EvaluationFlags::values |
                                 EvaluationFlags::gradients,
                               out);
      }
  };

  std::function<void(const MatrixFree<dim> &,
                     LinearAlgebra::distributed::Vector<double> &,
                     const LinearAlgebra::distributed::Vector<double> &,
                     const std::pair<unsigned int, unsigned int> &)>
    boundary_func =
      [](const MatrixFree<dim>                            &data,
         LinearAlgebra::distributed::Vector<double>       &out,
         const LinearAlgebra::distributed::Vector<double> &in,
         const std::pair<unsigned int, unsigned int>      &range) -> void {
    FEFaceEvaluation<dim, fe_degree> eval(data, true);

    for (unsigned int face = range.first; face < range.second; ++face)
      {
        eval.reinit(face);
        eval.gather_evaluate(in, EvaluationFlags::values);
        for (unsigned int q = 0; q < eval.n_q_points; ++q)
          {
            eval.submit_value(eval.quadrature_point(q).square() -
                                6. * eval.get_value(q),
                              q);
          }
        eval.integrate_scatter(EvaluationFlags::values, out);
      }
  };

  // compute reference result
  in.update_ghost_values();
  cell_func(matrix_free,
            ref,
            in,
            std::make_pair(0U, matrix_free.n_cell_batches()));
  boundary_func(matrix_free,
                ref,
                in,
                std::make_pair(matrix_free.n_inner_face_batches(),
                               matrix_free.n_inner_face_batches() +
                                 matrix_free.n_boundary_face_batches()));
  ref.compress(VectorOperation::add);
  in.zero_out_ghost_values();

  std::function<void(const MatrixFree<dim> &,
                     LinearAlgebra::distributed::Vector<double> &,
                     const LinearAlgebra::distributed::Vector<double> &,
                     const std::pair<unsigned int, unsigned int> &)>
    inner_face_func;


  // compute result through loop with various update settings. Note that we do
  // no compute on inner faces, so MatrixFree<dim>::DataAccessOnFaces::none
  // should be enough.
  matrix_free.loop(cell_func,
                   inner_face_func,
                   boundary_func,
                   test,
                   in,
                   true,
                   MatrixFree<dim>::DataAccessOnFaces::values,
                   MatrixFree<dim>::DataAccessOnFaces::values);

  deallog << "Number of dofs: " << dof_handler.n_dofs()
          << (overlap_communication ? " with overlap" : " without overlap")
          << std::endl;

  test -= ref;
  deallog << "Error loop 1: " << test.linfty_norm() << std::endl;
  matrix_free.loop(cell_func,
                   inner_face_func,
                   boundary_func,
                   test,
                   in,
                   true,
                   MatrixFree<dim>::DataAccessOnFaces::values,
                   MatrixFree<dim>::DataAccessOnFaces::values);
  test -= ref;
  deallog << "Error loop 2: " << test.linfty_norm() << std::endl;
  matrix_free.loop(cell_func,
                   inner_face_func,
                   boundary_func,
                   test,
                   in,
                   true,
                   MatrixFree<dim>::DataAccessOnFaces::none,
                   MatrixFree<dim>::DataAccessOnFaces::none);
  test -= ref;
  deallog << "Error loop 3: " << test.linfty_norm() << std::endl;
  matrix_free.loop(cell_func,
                   inner_face_func,
                   boundary_func,
                   test,
                   in,
                   true,
                   MatrixFree<dim>::DataAccessOnFaces::values,
                   MatrixFree<dim>::DataAccessOnFaces::none);
  test -= ref;
  deallog << "Error loop 4: " << test.linfty_norm() << std::endl;
  matrix_free.loop(cell_func,
                   inner_face_func,
                   boundary_func,
                   test,
                   in,
                   true,
                   MatrixFree<dim>::DataAccessOnFaces::unspecified,
                   MatrixFree<dim>::DataAccessOnFaces::unspecified);
  test -= ref;
  deallog << "Error loop 5: " << test.linfty_norm() << std::endl;

  // compute again, now only cell loop
  ref = 0;
  in.update_ghost_values();
  cell_func(matrix_free,
            ref,
            in,
            std::make_pair(0U, matrix_free.n_cell_batches()));
  ref.compress(VectorOperation::add);
  in.zero_out_ghost_values();

  matrix_free.cell_loop(cell_func, test, in, true);
  test -= ref;
  deallog << "Error cell loop: " << test.linfty_norm() << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  mpi_initlog();

  do_test<2, 1>(2, false);
  do_test<2, 1>(2, true);
  do_test<2, 1>(3, false);
  do_test<2, 1>(3, true);
  do_test<2, 2>(2, false);
  do_test<2, 2>(2, true);
  do_test<2, 2>(3, false);
  do_test<2, 2>(3, true);
  do_test<3, 2>(2, false);
  do_test<3, 2>(2, true);
}
