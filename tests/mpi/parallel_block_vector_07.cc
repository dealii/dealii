// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this BlockVector<Number>::multivector_inner_product_with_metric().
// Triangulation and Mass operator are the same as in
// matrix_free/mass_operator_01.cc

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/operators.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim, int fe_degree>
void
test(const unsigned int n_blocks = 5)
{
  using number = double;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  cell                                                   = tria.begin_active();
  for (; cell != endc; ++cell)
    if (cell->is_locally_owned())
      if (cell->center().norm() < 0.2)
        cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (dim < 3 && fe_degree < 2)
    tria.refine_global(2);
  else
    tria.refine_global(1);
  if (tria.begin(tria.n_levels() - 1)->is_locally_owned())
    tria.begin(tria.n_levels() - 1)->set_refine_flag();
  if (tria.last()->is_locally_owned())
    tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 10 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (cell->is_locally_owned())
          if (counter % (7 - i) == 0)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  const IndexSet &owned_set    = dof.locally_owned_dofs();
  const IndexSet  relevant_set = DoFTools::extract_locally_relevant_dofs(dof);

  AffineConstraints<double> constraints(owned_set, relevant_set);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  std::shared_ptr<MatrixFree<dim, number>> mf_data(
    new MatrixFree<dim, number>());
  {
    const QGauss<1>                                  quad(fe_degree + 2);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.tasks_block_size      = 7;
    mf_data->reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  MatrixFreeOperators::MassOperator<dim,
                                    fe_degree,
                                    fe_degree + 2,
                                    1,
                                    LinearAlgebra::distributed::Vector<number>>
    mf;
  mf.initialize(mf_data);
  mf.compute_diagonal();

  LinearAlgebra::distributed::BlockVector<number> left(n_blocks),
    right(n_blocks);
  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      mf_data->initialize_dof_vector(left.block(b));
      mf_data->initialize_dof_vector(right.block(b));
      left.block(b)  = 0.;
      right.block(b) = 0.;
      for (unsigned int i = 0; i < right.block(b).locally_owned_size(); ++i)
        {
          const unsigned int glob_index = owned_set.nth_index_in_set(i);
          if (constraints.is_constrained(glob_index))
            continue;
          right.block(b).local_element(i) = random_value<double>();
          left.block(b).local_element(i)  = random_value<double>();
        }
    }

  FullMatrix<number> metric(n_blocks, n_blocks);
  for (unsigned int i = 0; i < n_blocks; ++i)
    for (unsigned int j = 0; j < n_blocks; ++j)
      metric(i, j) = 1. + (i * 3 + j * 7);

  const double res = left.multivector_inner_product_with_metric(metric, right);

  double res2 = 0.;
  for (unsigned int i = 0; i < n_blocks; ++i)
    for (unsigned int j = 0; j < n_blocks; ++j)
      res2 += metric(i, j) * (left.block(i) * right.block(j));

  const double diff_norm = std::abs(res - res2);
  deallog << "Norm of difference: " << diff_norm << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test<2, 1>();
    }
  else
    {
      test<2, 1>();
    }
}
