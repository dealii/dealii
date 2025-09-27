// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests matrix-free partitioners for update_ghost_values and compress(add)


#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include "../tests.h"


template <int dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test(const std::vector<FiniteElement<dim, dim> *> &finite_elements)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  std::vector<std::shared_ptr<DoFHandler<dim, dim>>> dof_handlers(
    finite_elements.size());
  std::vector<const DoFHandler<dim, dim> *> dof_handlers_(dof_handlers.size());

  for (unsigned int i = 0; i < dof_handlers.size(); ++i)
    {
      dof_handlers[i] = std::make_shared<DoFHandler<dim, dim>>(tria);
      dof_handlers[i]->distribute_dofs(*finite_elements[i]);

      dof_handlers_[i] = dof_handlers[i].get();
    }

  MappingQ<dim> mapping(1);

  std::vector<Quadrature<1>> quads{QGauss<1>(finite_elements[0]->degree + 1)};

  AffineConstraints<Number> constraint;
  constraint.close();

  std::vector<const AffineConstraints<Number> *> constraints(
    dof_handlers.size());
  for (unsigned int i = 0; i < dof_handlers.size(); ++i)
    constraints[i] = &constraint;


  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData::none;
  additional_data.mapping_update_flags =
    update_gradients | update_JxW_values | update_quadrature_points;
  additional_data.mapping_update_flags_inner_faces =
    update_gradients | update_JxW_values | update_quadrature_points;
  additional_data.mapping_update_flags_boundary_faces =
    update_gradients | update_JxW_values | update_quadrature_points;
  additional_data.mapping_update_flags_faces_by_cells =
    update_gradients | update_JxW_values | update_quadrature_points;
  additional_data.hold_all_faces_to_owned_cells = true;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;

  matrix_free.reinit(
    mapping, dof_handlers_, constraints, quads, additional_data);

  deallog << "OK!" << std::endl;
}



int
main()
{
  initlog();

  {
    const unsigned int dim = 2;
    deallog.push("2d-dgq");
    FE_DGQ<dim> fe_dgq(1);
    test<dim>({&fe_dgq});
    deallog.pop();
  }

  {
    const unsigned int dim = 2;
    deallog.push("2d-q");
    FE_Q<dim> fe_q(1);
    test<dim>({&fe_q});
    deallog.pop();
  }

  {
    const unsigned int dim = 2;
    deallog.push("2d-dgq-q");
    FE_Q<dim>   fe_q(1);
    FE_DGQ<dim> fe_dgq(1);
    test<dim>({&fe_q, &fe_dgq});
    deallog.pop();
  }
}
