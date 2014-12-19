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




// DerivativeApproximation didn't work in parallel at all. This test verifies
// that it now does.


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
class Quadratic :
  public Function<dim>
{
public:
  double value (const Point<dim> &p, const unsigned int) const
    {
      return p*p;
    }
};


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

  // create a vector representing a function that is independent of the number
  // of processors involved
  TrilinosWrappers::MPI::Vector vec (dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  VectorTools::interpolate (dofh, Quadratic<dim>(), vec);
  vec.compress(VectorOperation::insert);

  // create such a vector with ghost elements and use it to compute
  // derivative information
  TrilinosWrappers::MPI::Vector vec_rel ( locally_relevant_set);
  vec_rel = vec;

  Vector<float> indicators(tr.n_active_cells());
  DerivativeApproximation::approximate_gradient  (dofh,
						  vec_rel,
						  indicators);

  // what we get must be a set of derivative indicators, one for each
  // cell of the distributed mesh. they need to be the same, no matter
  // how many processors we use. So, the sum of absolute values must
  // be the same for any processor number
  const double sum = Utilities::MPI::sum (indicators.l1_norm(), MPI_COMM_WORLD);
  if (myid == 0)
    deallog << sum << std::endl;
}



int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;

  test<2> ();
  test<3> ();
}

