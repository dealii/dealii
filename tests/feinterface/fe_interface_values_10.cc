// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check FEInterfaceViews::get_*_function_values/gradients...

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <fstream>
#include <type_traits>

#include "../tests.h"

template <int dim = 2>
void
test_jump_function()
{
  FE_DGP<dim> fe(0);

  Triangulation<dim> tria;
  DoFHandler<dim>    dof_handler(tria);
  Vector<double>     solution;

  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());

  std::vector<types::global_dof_index> local_dof_indices(1);
  // Populate with non-trivial values
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_indices(local_dof_indices);
      solution(local_dof_indices[0]) =
        static_cast<double>(cell->active_cell_index());
    }

  const QGauss<dim - 1>  face_quadrature(2);
  FEInterfaceValues<dim> fe_iv(fe,
                               face_quadrature,
                               update_values | update_gradients |
                                 update_quadrature_points |
                                 update_3rd_derivatives);

  std::vector<double> qp_jumps_global(face_quadrature.size());

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  for (; cell != dof_handler.end(); ++cell)
    for (const auto f : cell->face_indices())
      if (!cell->face(f)->at_boundary())
        {
          {
            fe_iv.reinit(cell,
                         f,
                         numbers::invalid_unsigned_int,
                         cell->neighbor(f),
                         cell->neighbor_of_neighbor(f),
                         numbers::invalid_unsigned_int);

            fe_iv.get_jump_in_function_values(solution, qp_jumps_global);

            double exact = cell->active_cell_index() * 1.0 -
                           cell->neighbor(f)->active_cell_index() * 1.0;

            for (unsigned int q = 0; q < face_quadrature.size(); ++q)
              Assert(std::fabs(qp_jumps_global[q] - exact) < 1e-15,
                     ExcNotImplemented());
          }
        }

  deallog << "OK" << std::endl;
}

template <int dim = 2>
void
test()
{
  FESystem<dim>                    fe(FE_DGP<dim>(4), 1);
  const FEValuesExtractors::Scalar extractor(0);

  Triangulation<dim> tria;
  DoFHandler<dim>    dof_handler(tria);
  Vector<double>     solution;

  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  dof_handler.distribute_dofs(fe);

  solution.reinit(dof_handler.n_dofs());

  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  // Populate with non-trivial values
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        solution(local_dof_indices[i]) =
          static_cast<double>(cell->active_cell_index() + 1);
    }

  const QGauss<dim - 1>  face_quadrature(2);
  FEInterfaceValues<dim> fe_iv(fe,
                               face_quadrature,
                               update_values | update_gradients |
                                 update_quadrature_points |
                                 update_3rd_derivatives);
  const unsigned int     n_face_q_points = face_quadrature.size();

  std::vector<double> jump_in_values(n_face_q_points);
  std::vector<double> exact_jump_in_values(n_face_q_points);

  std::vector<Tensor<1, dim>> jump_in_grads(n_face_q_points);
  std::vector<Tensor<1, dim>> exact_jump_in_grads(n_face_q_points);

  std::vector<Tensor<2, dim>> jump_in_hessians(n_face_q_points);
  std::vector<Tensor<2, dim>> exact_jump_in_hessians(n_face_q_points);

  std::vector<Tensor<3, dim>> jump_in_3rd_derivatives(n_face_q_points);
  std::vector<Tensor<3, dim>> exact_jump_in_3rd_derivatives(n_face_q_points);

  std::vector<double> average_of_values(n_face_q_points);
  std::vector<double> exact_average_of_values(n_face_q_points);

  std::vector<Tensor<1, dim>> average_of_grads(n_face_q_points);
  std::vector<Tensor<1, dim>> exact_average_of_grads(n_face_q_points);

  std::vector<Tensor<2, dim>> average_of_hessians(n_face_q_points);
  std::vector<Tensor<2, dim>> exact_average_of_hessians(n_face_q_points);

  std::vector<Tensor<3, dim>> average_of_3rd_derivatives(n_face_q_points);
  std::vector<Tensor<3, dim>> exact_average_of_3rd_derivatives(n_face_q_points);

  const double tol = 1e-15;

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  for (; cell != dof_handler.end(); ++cell)
    for (const auto f : cell->face_indices())
      if (!cell->face(f)->at_boundary())
        {
          {
            fe_iv.reinit(cell,
                         f,
                         numbers::invalid_unsigned_int,
                         cell->neighbor(f),
                         cell->neighbor_of_neighbor(f),
                         numbers::invalid_unsigned_int);

            fe_iv.get_jump_in_function_values(solution, jump_in_values);
            fe_iv[extractor].get_jump_in_function_values(solution,
                                                         exact_jump_in_values);

            fe_iv.get_jump_in_function_gradients(solution, jump_in_grads);
            fe_iv[extractor].get_jump_in_function_gradients(
              solution, exact_jump_in_grads);

            fe_iv.get_jump_in_function_hessians(solution, jump_in_hessians);
            fe_iv[extractor].get_jump_in_function_hessians(
              solution, exact_jump_in_hessians);

            fe_iv.get_jump_in_function_third_derivatives(
              solution, jump_in_3rd_derivatives);
            fe_iv[extractor].get_jump_in_function_third_derivatives(
              solution, exact_jump_in_3rd_derivatives);

            fe_iv.get_average_of_function_values(solution, average_of_values);
            fe_iv[extractor].get_average_of_function_values(
              solution, exact_average_of_values);

            fe_iv.get_average_of_function_gradients(solution, average_of_grads);
            fe_iv[extractor].get_average_of_function_gradients(
              solution, exact_average_of_grads);

            fe_iv.get_average_of_function_hessians(solution,
                                                   average_of_hessians);
            fe_iv[extractor].get_average_of_function_hessians(
              solution, exact_average_of_hessians);

            for (unsigned int q = 0; q < face_quadrature.size(); ++q)
              {
                Assert(std::fabs(jump_in_values[q] - exact_jump_in_values[q]) <
                         tol,
                       ExcNotImplemented());
                Assert((jump_in_grads[q] - exact_jump_in_grads[q]).norm() < tol,
                       ExcNotImplemented());
                Assert((jump_in_hessians[q] - exact_jump_in_hessians[q])
                           .norm() < tol,
                       ExcNotImplemented());
                Assert((jump_in_3rd_derivatives[q] -
                        exact_jump_in_3rd_derivatives[q])
                           .norm() < tol,
                       ExcNotImplemented());

                Assert(std::fabs(average_of_values[q] -
                                 exact_average_of_values[q]) < tol,
                       ExcNotImplemented());
                Assert((average_of_grads[q] - exact_average_of_grads[q])
                           .norm() < tol,
                       ExcNotImplemented());
                Assert((average_of_hessians[q] - exact_average_of_hessians[q])
                           .norm() < tol,
                       ExcNotImplemented());
              }
          }
        }

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.push("Test jump in function values");
  {
    test_jump_function();
  }
  deallog.pop();
  deallog << "OK" << std::endl;

  deallog.push("Test all functions");
  {
    test();
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
