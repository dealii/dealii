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



// check ConstraintMatrix for a distributed mesh on hyper shell with both
// hanging nodes and normal-flux constraints on the outer boundary of the
// shell.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <sstream>

template<int dim>
void test()
{
  Assert (dim == 3, ExcNotImplemented());
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);

  // create hypershell mesh and refine some cells. use some large numbers
  // which make round-off influence more pronounced
  const double R0      = 6371000.-2890000.;     /* m          */
  const double R1      = 6371000.-  35000.;     /* m          */
  GridGenerator::hyper_shell (triangulation,
                              Point<dim>(),
                              R0,
                              R1,
                              96, true);
  triangulation.refine_global (2);

  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (!cell->is_ghost() && !cell->is_artificial())
      if (cell->center()[2] > 0.75 * R1)
        {
          cell->set_refine_flag();
          for (unsigned int c=0; c<8; ++c)
            cell->parent()->child(c)->set_refine_flag();
        }

  {
    unsigned int n_flagged_cells = 0;
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      if (!cell->is_ghost() && !cell->is_artificial())
        if (cell->refine_flag_set())
          ++n_flagged_cells;

    unsigned int global_f_c = 0;
    MPI_Allreduce (&n_flagged_cells, &global_f_c, 1, MPI_UNSIGNED,
                   MPI_SUM, MPI_COMM_WORLD);
    if (myid == 0)
      deallog << "# flagged cells = " << global_f_c << std::endl;
  }

  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement ();

  if (myid == 0)
    deallog <<  "#cells = " << triangulation.n_global_active_cells()
            << std::endl;

  // create FE_System and fill in no-normal flux
  // conditions on boundary 1 (outer)
  static const FESystem<dim> fe(FE_Q<dim> (1), dim);
  DoFHandler<dim> dofh(triangulation);
  dofh.distribute_dofs (fe);
  if (myid == 0)
    deallog <<  "#dofs = " << dofh.locally_owned_dofs().size()
            << std::endl;

  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);

  ConstraintMatrix constraints;
  constraints.reinit(relevant_set);
  DoFTools::make_hanging_node_constraints (dofh, constraints);
  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (1);
  const unsigned int degree = 1;
  VectorTools::compute_normal_flux_constraints (dofh, 0,
                                                no_normal_flux_boundaries,
                                                constraints,
                                                MappingQ<dim>(degree));
  constraints.close();

  // print the number of constraints. since
  // processors might write info in different
  // orders, copy all numbers to root processor
  std::vector<unsigned int> n_constraints_glob (numprocs);
  unsigned int n_constraints = constraints.n_constraints();
  MPI_Gather (&n_constraints, 1, MPI_UNSIGNED,
              &n_constraints_glob[0], 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if (myid == 0)
    for (unsigned int i=0; i<numprocs; ++i)
      deallog <<  "#constraints on " << i << ": " << n_constraints_glob[i]
              << std::endl;

  // dummy assembly: put 1 in all components of
  // the vector
  TrilinosWrappers::MPI::Vector vector;
  vector.reinit (dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  {
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    Vector<double> local_vector (dofs_per_cell);
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      local_vector(i) = 1.;
    typename DoFHandler<dim>::active_cell_iterator
    cell = dofh.begin_active(),
    endc = dofh.end();
    for (; cell!=endc; ++cell)
      if (cell->subdomain_id() == triangulation.locally_owned_subdomain())
        {
          cell->get_dof_indices (local_dof_indices);
          constraints.distribute_local_to_global (local_vector,
                                                  local_dof_indices,
                                                  vector);
        }
    vector.compress (VectorOperation::add);
  }

  // now check that no entries were generated
  // for constrained entries on the locally
  // owned range.
  const std::pair<unsigned int,unsigned int> range = vector.local_range();
  for (unsigned int i=range.first; i<range.second; ++i)
    if (constraints.is_constrained(i))
      Assert (vector(i)==0, ExcInternalError());

  if (myid==0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  {
    Utilities::MPI::MPI_InitFinalize mpi_init (argc, argv);
    unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


    deallog.push(Utilities::int_to_string(myid));

    if (myid == 0)
      {
        std::ofstream logfile("output");
        deallog.attach(logfile);
        deallog.depth_console(0);
        deallog.threshold_double(1.e-10);

        deallog.push("3d");
        test<3>();
        deallog.pop();
      }
    else
      test<3>();
  }

  return 0;
}
