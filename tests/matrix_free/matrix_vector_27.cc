// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// same as matrix_vector_26.cc but for level operator

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

#include "matrix_vector_mf.h"



template <int dim, int fe_degree>
void
test()
{
  using number = double;

  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  unsigned int max_level = tria.n_global_levels() - 1;

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs();

  IndexSet owned_set = dof.locally_owned_dofs();
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs(dof, relevant_set);

  AffineConstraints<double> constraints(relevant_set);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  // std::cout << "Number of cells: " << tria.n_global_active_cells() <<
  // std::endl; std::cout << "Number of degrees of freedom: " << dof.n_dofs() <<
  // std::endl; std::cout << "Number of constraints: " <<
  // constraints.n_constraints() << std::endl;

  MatrixFree<dim, number> mf_data;
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.tasks_block_size      = 7;
    data.store_ghost_cells     = true;
    data.mg_level              = max_level;
    mf_data.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  }

  // Check if every ghosted cell is included in MatrixFree storage
  {
    using IndexLevel = std::pair<unsigned int, unsigned int>;
    std::set<std::pair<unsigned int, unsigned int>> ghost_cells;



    for (unsigned int batch = 0;
         batch < mf_data.n_cell_batches() + mf_data.n_ghost_cell_batches();
         ++batch)
      {
        for (unsigned int lane = 0;
             lane < mf_data.n_active_entries_per_cell_batch(batch);
             ++lane)
          {
            const typename Triangulation<dim>::level_cell_iterator cell =
              mf_data.get_cell_iterator(batch, lane);
            IndexLevel index_level(cell->index(), cell->level());
            ghost_cells.insert(index_level);
          }
      }

    for (const auto &cell : tria.cell_iterators_on_level(max_level))
      if (cell->is_ghost())
        {
          IndexLevel index_level(cell->index(), cell->level());

          AssertThrow(ghost_cells.find(index_level) != ghost_cells.end(),
                      ExcMessage(
                        "MatrixFree does not store all ghost cells. "
                        "Index: " +
                        std::to_string(index_level.first) +
                        ", level: " + std::to_string(index_level.second)));
        }
  }

  MatrixFreeTest<dim,
                 fe_degree,
                 number,
                 LinearAlgebra::distributed::Vector<number>>
                                             mf(mf_data);
  LinearAlgebra::distributed::Vector<number> in, out, ref;
  mf_data.initialize_dof_vector(in);
  out.reinit(in);
  ref.reinit(in);

  for (unsigned int i = 0; i < in.locally_owned_size(); ++i)
    {
      const unsigned int glob_index = owned_set.nth_index_in_set(i);
      if (constraints.is_constrained(glob_index))
        continue;
      in.local_element(i) = random_value<double>();
    }

  mf.vmult(out, in);


  // assemble trilinos sparse matrix with
  // (\nabla v, \nabla u) + (v, 10 * u) for
  // reference
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
      if (cell->is_locally_owned())
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

  sparse_matrix.vmult(ref, in);
  out -= ref;
  const double diff_norm = out.linfty_norm();

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
      test<2, 2>();
      deallog.pop();

      deallog.push("3d");
      test<3, 1>();
      test<3, 2>();
      deallog.pop();
    }
  else
    {
      test<2, 1>();
      test<2, 2>();
      test<3, 1>();
      test<3, 2>();
    }
}
