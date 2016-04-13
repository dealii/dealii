// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2015 by the deal.II authors
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


// Document bug in determining level subdomain ids in certain cases. This used
// to cause deadlocks because the list of ghost neighbors was not symmetric.
//
// problematic cases: 2d: 6 and 13 procs, 3d 20 procs (original bug report from Martin)

#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
void print(parallel::distributed::Triangulation<dim> &tr)
{
  deallog << "*****" << std::endl;
  for (typename parallel::distributed::Triangulation<dim>::cell_iterator cell = tr.begin();
       cell != tr.end(); ++cell)
    {
      if (cell->level_subdomain_id() != numbers::artificial_subdomain_id)
        deallog << "cell=" << cell->id()
                << " level_subdomain_id=" << cell->level_subdomain_id()
                << std::endl;
    }
}


template <int dim>
void do_test ()
{
  FE_Q<dim> fe(1);
  deallog << "Testing " << fe.get_name() << std::endl << std::endl;
  parallel::distributed::Triangulation<dim> triangulation
  (MPI_COMM_WORLD,
   Triangulation<dim>::none,
   parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  GridGenerator::subdivided_hyper_cube (triangulation, 1, -1, 1);
  triangulation.refine_global(2);
  for (typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (cell->is_locally_owned() &&
        cell->center().norm() < 0.55)
      cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  print(triangulation);

  for (typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (cell->is_locally_owned() &&
        cell->center().norm() > 0.3 && cell->center().norm() < 0.42)
      cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  print(triangulation);

  for (typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (cell->is_locally_owned() &&
        cell->center().norm() > 0.335 && cell->center().norm() < 0.39)
      cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  print(triangulation);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs(fe);
}


int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll log;

  do_test<2>();
  do_test<3>();
  return 0;
}
