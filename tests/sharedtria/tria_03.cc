// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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
//

// p::s::Tria n_active_cells() and n_global_active_cells() with artificial cells.

#include "../tests.h"
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include <set>
#include <cstdio>


template<int dim>
bool
pred_d(const typename parallel::shared::Triangulation<dim>::active_cell_iterator &cell)
{
  return (cell->center()(0) < 0.5 &&
          cell->center()(1) < 0.5);
}



template <int dim>
void test (const unsigned int flag)
{
  // Setup system
  parallel::shared::Triangulation<dim> triangulation(MPI_COMM_WORLD, dealii::Triangulation<dim>::none, true);

  GridGenerator::hyper_rectangle (triangulation,
                                  Point<dim>(0,0),
                                  Point<dim>(1,1));

  triangulation.refine_global(2);

  // Extra refinement to generate hanging nodes
  for (typename parallel::distributed::Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (cell->is_locally_owned() &&
        ( (flag==1 &&  pred_d<dim>(cell)) ||
          (flag==2 && !pred_d<dim>(cell))  )
       )
      cell->set_refine_flag ();

  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement ();

  deallog << "n_cells=" << triangulation.n_cells() << std::endl
          << "n_active_cells="<<triangulation.n_active_cells() << std::endl
          << "n_global_active_cells="<<triangulation.n_global_active_cells() << std::endl;

}

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  MPILogInitAll log;

  deallog.push(Utilities::int_to_string(myid));

  //test<2>(0);
  test<2>(1);
  //test<2>(2);

  deallog.pop();
}
