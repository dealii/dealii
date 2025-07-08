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

// Same as matrix_free_matrix_vector_25 test but we only perform the operator
// evaluation on the cells that have a fe_index equal to 0

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include "deal.II/grid/filtered_iterator.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

#include "matrix_vector_device_mf.h"



template <int dim, int fe_degree>
void
test()
{
  using Number = double;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
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

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim>(fe_degree));
  fe_collection.push_back(FE_Nothing<dim>());

  DoFHandler<dim> dof(tria);
  // Set FE_Index
  for (auto &cell_ : dof.active_cell_iterators())
    if (cell_->is_locally_owned())
      {
        if (cell_->center()[0] < 0.5)
          cell_->set_active_fe_index(0);
        else
          cell_->set_active_fe_index(1);
      }
  dof.distribute_dofs(fe_collection);

  const IndexSet &owned_set    = dof.locally_owned_dofs();
  const IndexSet  relevant_set = DoFTools::extract_locally_relevant_dofs(dof);

  AffineConstraints<double> constraints(owned_set, relevant_set);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  MappingQ<dim>                     mapping(fe_degree);
  Portable::MatrixFree<dim, Number> mf_data;
  const QGauss<1>                   quad(fe_degree + 1);
  typename Portable::MatrixFree<dim, Number>::AdditionalData additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients |
                                         update_JxW_values |
                                         update_quadrature_points;
  IteratorFilters::ActiveFEIndexEqualTo filter(0, true);
  mf_data.reinit(mapping, dof, constraints, quad, filter, additional_data);

  const unsigned int coef_size =
    tria.n_locally_owned_active_cells() * std::pow(fe_degree + 1, dim);
  MatrixFreeTest<
    dim,
    fe_degree,
    Number,
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>>
    mf(mf_data, coef_size);
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> in_dev;
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Default> out_dev;
  mf_data.initialize_dof_vector(in_dev);
  mf_data.initialize_dof_vector(out_dev);

  LinearAlgebra::ReadWriteVector<Number> rw_in(owned_set);
  for (unsigned int i = 0; i < in_dev.locally_owned_size(); ++i)
    {
      const unsigned int glob_index = owned_set.nth_index_in_set(i);
      if (constraints.is_constrained(glob_index))
        continue;
      rw_in.local_element(i) = random_value<double>();
    }
  in_dev.import_elements(rw_in, VectorOperation::insert);
  mf.vmult(out_dev, in_dev);

  LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> out_host(
    owned_set, MPI_COMM_WORLD);
  LinearAlgebra::ReadWriteVector<Number> rw_out(owned_set);
  rw_out.import_elements(out_dev, VectorOperation::insert);
  out_host.import_elements(rw_out, VectorOperation::insert);

  // assemble trilinos sparse matrix with
  // (\nabla v, \nabla u) + (v, 10 * u) for
  // reference
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> in_host(
    owned_set, MPI_COMM_WORLD);
  in_host.import_elements(rw_in, VectorOperation::insert);
  LinearAlgebra::distributed::Vector<Number, MemorySpace::Host> ref(
    owned_set, MPI_COMM_WORLD);
  TrilinosWrappers::SparseMatrix sparse_matrix;
  {
    TrilinosWrappers::SparsityPattern csp(owned_set, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof,
                                    csp,
                                    constraints,
                                    true,
                                    Utilities::MPI::this_mpi_process(
                                      MPI_COMM_WORLD));
    csp.compress();
    sparse_matrix.reinit(csp);
  }
  {
    QGauss<dim> quadrature_formula(fe_degree + 1);

    FEValues<dim> fe_values(dof.get_fe(),
                            quadrature_formula,
                            update_values | update_gradients |
                              update_JxW_values);

    const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(),
                                                   endc = dof.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned() && cell->active_fe_index() == 0)
        {
          cell_matrix = 0;
          fe_values.reinit(cell);

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) +=
                    ((fe_values.shape_grad(i, q_point) *
                        fe_values.shape_grad(j, q_point) +
                      10. * fe_values.shape_value(i, q_point) *
                        fe_values.shape_value(j, q_point)) *
                     fe_values.JxW(q_point));
              }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 sparse_matrix);
        }
  }
  sparse_matrix.compress(VectorOperation::add);

  sparse_matrix.vmult(ref, in_host);
  out_host -= ref;

  const double diff_norm = out_host.linfty_norm();

  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
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

      deallog.push("2d");
      test<2, 1>();
      deallog.pop();

      deallog.push("3d");
      test<3, 1>();
      test<3, 2>();
      deallog.pop();
    }
  else
    {
      test<2, 1>();
      test<3, 1>();
      test<3, 2>();
    }
}
