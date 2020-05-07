// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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



// check FEFaceEvaluation::read_dof_values_plain on a continuous finite
// element

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim, int fe_degree>
void
do_test(const unsigned int n_refine)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0.1, 0.96);
  tria.refine_global(n_refine);

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  MatrixFree<dim>                          matrix_free;
  typename MatrixFree<dim>::AdditionalData add_data;
  add_data.mapping_update_flags = update_values | update_gradients |
                                  update_JxW_values | update_quadrature_points;
  add_data.mapping_update_flags_boundary_faces =
    update_values | update_JxW_values | update_quadrature_points;
  matrix_free.reinit(dof_handler,
                     constraints,
                     QGauss<1>(fe.degree + 1),
                     add_data);

  LinearAlgebra::distributed::Vector<double> in, test;
  matrix_free.initialize_dof_vector(in);
  for (unsigned int i = 0; i < in.get_partitioner()->local_size(); ++i)
    in.local_element(i) = in.get_partitioner()->local_to_global(i);
  matrix_free.initialize_dof_vector(test);

  std::function<void(const MatrixFree<dim> &,
                     LinearAlgebra::distributed::Vector<double> &,
                     const LinearAlgebra::distributed::Vector<double> &,
                     const std::pair<unsigned int, unsigned int> &)>
    cell_func = [](const MatrixFree<dim> &                           data,
                   LinearAlgebra::distributed::Vector<double> &      out,
                   const LinearAlgebra::distributed::Vector<double> &in,
                   const std::pair<unsigned int, unsigned int> &range) -> void {
    FEEvaluation<dim, fe_degree> eval(data);

    for (unsigned int cell = range.first; cell < range.second; ++cell)
      {
        eval.reinit(cell);
        eval.read_dof_values_plain(in);
        eval.evaluate(false, true);
        for (unsigned int q = 0; q < eval.n_q_points; ++q)
          {
            eval.submit_gradient(eval.get_gradient(q), q);
            eval.submit_value(eval.quadrature_point(q).square(), q);
          }
        eval.integrate_scatter(true, true, out);
      }
  };

  std::function<void(const MatrixFree<dim> &,
                     LinearAlgebra::distributed::Vector<double> &,
                     const LinearAlgebra::distributed::Vector<double> &,
                     const std::pair<unsigned int, unsigned int> &)>
    boundary_func =
      [](const MatrixFree<dim> &                           data,
         LinearAlgebra::distributed::Vector<double> &      out,
         const LinearAlgebra::distributed::Vector<double> &in,
         const std::pair<unsigned int, unsigned int> &     range) -> void {
    FEFaceEvaluation<dim, fe_degree> eval(data, true);

    for (unsigned int face = range.first; face < range.second; ++face)
      {
        eval.reinit(face);
        eval.read_dof_values_plain(in);
        eval.evaluate(true, false);
        for (unsigned int q = 0; q < eval.n_q_points; ++q)
          {
            eval.submit_value(eval.quadrature_point(q).square() -
                                6. * eval.get_value(q),
                              q);
          }
        eval.integrate_scatter(true, false, out);
      }
  };

  std::function<void(const MatrixFree<dim> &,
                     LinearAlgebra::distributed::Vector<double> &,
                     const LinearAlgebra::distributed::Vector<double> &,
                     const std::pair<unsigned int, unsigned int> &)>
    inner_face_func;

  matrix_free.loop(cell_func,
                   inner_face_func,
                   boundary_func,
                   test,
                   in,
                   true,
                   MatrixFree<dim>::DataAccessOnFaces::values,
                   MatrixFree<dim>::DataAccessOnFaces::values);

  deallog << "L2 norm of result: " << test.l2_norm() << std::endl;
}



int
main()
{
  initlog();

  do_test<2, 1>(0);
  do_test<2, 1>(1);
  do_test<2, 1>(2);
  do_test<3, 1>(0);
  do_test<3, 1>(1);
  do_test<3, 1>(2);
  do_test<3, 2>(1);
}
