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

// test FEValuesBase::get_function_values()
// when used with complex vector;


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

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include <deal.II/lac/parallel_vector.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

using namespace dealii;

static const unsigned int dim = 2;

void test ()
{

  MPI_Comm mpi_communicator (MPI_COMM_WORLD);
  const unsigned int n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout (std::cout,
                            (Utilities::MPI::this_mpi_process(mpi_communicator)
                             == 0));

  parallel::distributed::Triangulation<dim> tria(mpi_communicator,
                                                 typename Triangulation<dim>::MeshSmoothing
                                                 (Triangulation<dim>::smoothing_on_refinement |
                                                  Triangulation<dim>::smoothing_on_coarsening));


  GridGenerator::hyper_cube (tria, -1,0);
  tria.refine_global (2);


  const unsigned int poly_degree = 1;
  FE_Q<dim> fe(poly_degree);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs ();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dof_handler,locally_relevant_dofs);

  PETScWrappers::MPI::Vector vector, vector_locally_relevant;
  parallel::distributed::Vector< double > vector_Re, vector_Re_locally_relevant,
           vector_Im, vector_Im_locally_relevant;
  vector.reinit(locally_owned_dofs, mpi_communicator);
  vector_locally_relevant.reinit (locally_owned_dofs,
                                  locally_relevant_dofs,mpi_communicator);

  vector_Re.reinit(locally_owned_dofs, mpi_communicator);
  vector_Re_locally_relevant.reinit (locally_owned_dofs,
                                     locally_relevant_dofs,mpi_communicator);

  vector_Im.reinit(locally_owned_dofs, mpi_communicator);
  vector_Im_locally_relevant.reinit (locally_owned_dofs,
                                     locally_relevant_dofs,mpi_communicator);

  const types::global_dof_index n_local_dofs = locally_owned_dofs.n_elements();

  ConstraintMatrix constraints;
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  constraints.close ();
  //set vector:
  unsigned int myid = Utilities::MPI::this_mpi_process (mpi_communicator);
  for (unsigned int i = 0; i < locally_owned_dofs.n_elements(); i++)
    {
      const PetscScalar val = 1.0+myid+ (myid+i%2)*2.0*PETSC_i;
      vector   ( locally_owned_dofs.nth_index_in_set(i)) = val;
      vector_Re( locally_owned_dofs.nth_index_in_set(i)) = PetscRealPart (val);
      vector_Im( locally_owned_dofs.nth_index_in_set(i)) = PetscImaginaryPart (val);

    }
  vector.compress(VectorOperation::insert);
  vector_Re.compress(VectorOperation::insert);
  vector_Im.compress(VectorOperation::insert);

  vector_locally_relevant = vector;
  vector_Re_locally_relevant = vector_Re;
  vector_Im_locally_relevant = vector_Im;

  //test get_function_values:
  {
    //a separate quadrature formula
    //enough for mass and kinetic matrices assembly
    QGauss<dim> quadrature_formula(poly_degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    std::vector<PetscScalar> values(n_q_points);
    std::vector<double> values_Re(n_q_points),
        values_Im(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active (),
    endc = dof_handler.end ();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values.get_function_values (vector_locally_relevant,    values);
          fe_values.get_function_values (vector_Re_locally_relevant, values_Re);
          fe_values.get_function_values (vector_Im_locally_relevant, values_Im);

          for (unsigned int q=0; q<n_q_points; ++q)
            {
              AssertThrow ( PetscRealPart(values[q])      == values_Re[q], ExcInternalError() );
              AssertThrow ( PetscImaginaryPart(values[q]) == values_Im[q], ExcInternalError() );
            }
        }
  }

  deallog << "OK" << std::endl;
}


int main (int argc, char *argv[])
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  {
    test ();
  }
}
