// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2016 by the deal.II authors
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

// test that matrix_scalar_product of a symmetric matrix
// applied to the same vector result in a real number


#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

using namespace dealii;

template<int dim>
void test (const unsigned int poly_degree = 1)
{

  MPI_Comm mpi_communicator (MPI_COMM_WORLD);
  const unsigned int n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator));

  parallel::distributed::Triangulation<dim> tria(mpi_communicator,
                                                 typename Triangulation<dim>::MeshSmoothing
                                                 (Triangulation<dim>::smoothing_on_refinement |
                                                  Triangulation<dim>::smoothing_on_coarsening));

  GridGenerator::hyper_cube (tria, -1,0);
  tria.refine_global (3);

  FE_Q<dim> fe(poly_degree);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs ();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dof_handler,locally_relevant_dofs);

  PETScWrappers::MPI::Vector vector;
  PETScWrappers::MPI::SparseMatrix mass_matrix;


  vector.reinit(locally_owned_dofs, mpi_communicator);

  ConstraintMatrix constraints;
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  constraints.close ();

  DynamicSparsityPattern dsp (locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dof_handler, dsp,
                                   constraints, false);
  SparsityTools::distribute_sparsity_pattern (dsp,
                                              dof_handler.n_locally_owned_dofs_per_processor(),
                                              mpi_communicator,
                                              locally_relevant_dofs);

  mass_matrix.reinit (locally_owned_dofs,
                      locally_owned_dofs,
                      dsp,
                      mpi_communicator);

  //assemble mass matrix:
  mass_matrix = PetscScalar();
  {
    QGauss<dim> quadrature_formula(poly_degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<PetscScalar> cell_mass_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active (),
    endc = dof_handler.end ();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          cell_mass_matrix    = PetscScalar();

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                {
                  cell_mass_matrix (i, j)
                  += (fe_values.shape_value (i, q_point) *
                      fe_values.shape_value (j, q_point)
                     ) * fe_values.JxW (q_point);
                }

          cell->get_dof_indices (local_dof_indices);
          constraints.distribute_local_to_global (cell_mass_matrix,
                                                  local_dof_indices,
                                                  mass_matrix);
        }

    mass_matrix.compress (VectorOperation::add);
  }

  for (unsigned int i = 0; i < locally_owned_dofs.n_elements(); i++)
    {
      double re = 0,
             im = 0;
      if ( i % 2 )
        {
          re = 1.0*i;
          im = 1.0*(this_mpi_process+1);
        }
      else if (i % 3)
        {
          re = 0.0;
          im = -1.0*(this_mpi_process+1);
        }
      else
        {
          re = 3.0*i;
          im = 0.0;
        }
      const PetscScalar val = re+im*PETSC_i;
      vector   ( locally_owned_dofs.nth_index_in_set(i)) = val;
    }
  vector.compress(VectorOperation::insert);
  constraints.distribute(vector);

  PETScWrappers::MPI::Vector tmp(vector);
  mass_matrix.vmult (tmp, vector);
  //mass_matrix.Tvmult (tmp, vector);

  const std::complex<double> norm1 = vector*tmp;
  deallog<<"(conj(vector),M vector): "<<std::endl;
  deallog<<"real part:      "<<PetscRealPart(norm1)<<std::endl;
  deallog<<"imaginary part: "<<PetscImaginaryPart(norm1)<<std::endl;

  const std::complex<double> norm2 =
    mass_matrix.matrix_scalar_product(vector, vector);
  deallog<<"matrix_scalar_product(vec,vec): "<<std::endl;
  deallog<<"real part:      "<<PetscRealPart(norm2)<<std::endl;
  deallog<<"imaginary part: "<<PetscImaginaryPart(norm2)<<std::endl;

  const std::complex<double> norm3 =
    mass_matrix.matrix_norm_square(vector);
  deallog<<"matrix_norm_square(vec): "<<std::endl;
  deallog<<"real part:      "<<PetscRealPart(norm3)<<std::endl;
  deallog<<"imaginary part: "<<PetscImaginaryPart(norm3)<<std::endl;
}


int main (int argc, char *argv[])
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
      test<2>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
