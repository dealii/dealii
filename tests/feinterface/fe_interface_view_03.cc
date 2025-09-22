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
test()
{
  FESystem<dim>                    fe(FE_DGP<dim>(4), FE_DGP<dim>(2) ^ 2);
  const FEValuesExtractors::Scalar scalar_extractor(0);
  const FEValuesExtractors::Vector vector_extractor(1);

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
                               update_gradients | update_quadrature_points);

  const double tol = 1e-15;

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();

  deallog.push("Test scalar values");
  {
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

              for (unsigned int q = 0; q < face_quadrature.size(); ++q)
                {
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      const Tensor<1, dim> gradient_here =
                        fe_iv[scalar_extractor].gradient(true, i, q);
                      const Tensor<1, dim> gradient_there =
                        fe_iv[scalar_extractor].gradient(false, i, q);

                      const Tensor<1, dim> manual_jump_in_grads =
                        gradient_here - gradient_there;
                      const Tensor<1, dim> given_jump_in_grads =
                        fe_iv[scalar_extractor].jump_in_gradients(i, q);
                      Assert((manual_jump_in_grads - given_jump_in_grads)
                                 .norm() < tol,
                             ExcNotImplemented());

                      const Tensor<1, dim> manual_average_of_grads =
                        0.5 * (gradient_here + gradient_there);
                      const Tensor<1, dim> given_average_of_grads =
                        fe_iv[scalar_extractor].average_of_gradients(i, q);
                      Assert((manual_average_of_grads - given_average_of_grads)
                                 .norm() < tol,
                             ExcNotImplemented());
                    }
                }
            }
          }
  }
  deallog << "OK" << std::endl;
  deallog.pop();

  deallog.push("Test vector values");
  {
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

              for (unsigned int q = 0; q < face_quadrature.size(); ++q)
                {
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      const Tensor<2, dim> gradient_here =
                        fe_iv[vector_extractor].gradient(true, i, q);
                      const Tensor<2, dim> gradient_there =
                        fe_iv[vector_extractor].gradient(false, i, q);

                      const Tensor<2, dim> manual_jump_in_grads =
                        gradient_here - gradient_there;
                      const Tensor<2, dim> given_jump_in_grads =
                        fe_iv[vector_extractor].jump_in_gradients(i, q);
                      Assert((manual_jump_in_grads - given_jump_in_grads)
                                 .norm() < tol,
                             ExcNotImplemented());

                      const Tensor<2, dim> manual_average_of_grads =
                        0.5 * (gradient_here + gradient_there);
                      const Tensor<2, dim> given_average_of_grads =
                        fe_iv[vector_extractor].average_of_gradients(i, q);
                      Assert((manual_average_of_grads - given_average_of_grads)
                                 .norm() < tol,
                             ExcNotImplemented());
                    }
                }
            }
          }
  }
  deallog << "OK" << std::endl;
  deallog.pop();
}

int
main(int argc, char **argv)
{
  initlog();
  test();
}
