// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// document bug in the DerivateApproximation in parallel
/*
 --------------------------------------------------------
 An error occurred in line <3349> of file </ssd/deal-trunk/deal.II/include/deal.II/dofs/dof_accessor.templates.h> in function
     void dealii::DoFCellAccessor<DH, lda>::get_dof_values(const InputVector&, ForwardIterator, ForwardIterator) const [with InputVector = dealii::TrilinosWrappers::MPI::Vector, ForwardIterator = double*, DH = dealii::DoFHandler<2>, bool level_dof_access = false]
 The violated condition was:
     this->is_artificial() == false
 The name and call sequence of the exception was:
     ExcMessage ("Can't ask for DoF indices on artificial cells.")
 Additional Information:
 Can't ask for DoF indices on artificial cells.
 --------------------------------------------------------

 */

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/derivative_approximation.h>

#include <fstream>
#include <sstream>



template <int dim>
void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int n_processes = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  IndexSet locally_relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh,
                                           locally_relevant_set);

  TrilinosWrappers::MPI::Vector vec (dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  for (unsigned int i=vec.local_range().first; i<vec.local_range().second; ++i)
    vec(i) = i;
  vec.compress(VectorOperation::insert);

  TrilinosWrappers::MPI::Vector vec_rel ( locally_relevant_set);
  vec_rel = vec;

  MappingQ1<dim> mapping;
  Vector<float> indicators(tr.n_active_cells());
  DerivativeApproximation::approximate_gradient  (mapping,
						  dofh,
						  vec_rel,
						  indicators);

  // we got here, so no exception.
  if (myid == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;

  test<2>();
}
