// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// this tests the correctness of distribute_local_to_globa() for DiagonalMatrix

#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>




template <int dim, int fe_degree>
void test (const bool hanging_nodes = true)
{
  typedef double number;

  parallel::distributed::Triangulation<dim> tria (MPI_COMM_WORLD);
  GridGenerator::hyper_cube (tria);
  tria.refine_global(1);

  if (hanging_nodes)
    {
      if (tria.begin_active()->is_locally_owned())
        tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement ();
    }

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs(fe);

  IndexSet owned_set = dof.locally_owned_dofs();
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dof, relevant_set);

  ConstraintMatrix constraints (relevant_set);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values (dof, 0, ZeroFunction<dim>(),
                                            constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;
  //std::cout << "Number of cells: " << tria.n_global_active_cells() << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
  //std::cout << "Number of constraints: " << constraints.n_constraints() << std::endl;


  // assemble trilinos sparse matrix
  TrilinosWrappers::SparseMatrix sparse_matrix;
  {
    TrilinosWrappers::SparsityPattern csp (owned_set, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern (dof, csp, constraints, true,
                                     Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
    csp.compress();
    sparse_matrix.reinit (csp);
  }

  DiagonalMatrix<LinearAlgebra::distributed::Vector<double> > diagonal_matrix;
  {
    LinearAlgebra::distributed::Vector<double> dummy;
    dummy.reinit(owned_set, relevant_set, MPI_COMM_WORLD);
    diagonal_matrix.reinit(dummy);
  }

  {
    QGauss<dim>  quadrature_formula(fe_degree+1);

    FEValues<dim> fe_values (dof.get_fe(), quadrature_formula,
                             update_gradients | update_JxW_values);

    const unsigned int   dofs_per_cell = dof.get_fe().dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof.begin_active(),
    endc = dof.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          fe_values.reinit (cell);

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                       fe_values.shape_grad(j,q_point)) *
                                      fe_values.JxW(q_point);
              }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global (cell_matrix,
                                                  local_dof_indices,
                                                  sparse_matrix);
          constraints.distribute_local_to_global (cell_matrix,
                                                  local_dof_indices,
                                                  diagonal_matrix);
        }
  }
  sparse_matrix.compress(VectorOperation::add);
  diagonal_matrix.compress(VectorOperation::add);

  // compare elements:
  for (unsigned int i=0; i<owned_set.n_elements(); ++i)
    {
      const unsigned int glob_index =
        owned_set.nth_index_in_set (i);
      if (constraints.is_constrained(glob_index))
        continue;

      const double d = diagonal_matrix(glob_index, glob_index);
      const double t = sparse_matrix.diag_element(glob_index);
      if ( std::abs(d - t)/t > 1e-10 )
        {
          deallog << glob_index << " " << d << " != " << t << std::endl;
        }
    }

  deallog << "Ok." << std::endl;
}


int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  mpi_initlog();

  deallog.push("2d");
  test<2,1>();
  test<2,2>();
  deallog.pop();

  deallog.push("3d");
  test<3,1>();
  test<3,2>();
  deallog.pop();
}

