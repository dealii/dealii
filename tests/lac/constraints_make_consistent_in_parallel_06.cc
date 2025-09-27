// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test AffineConstraints::make_consistent_in_parallel() for the case where
// constraints need to be combined between different ranks. We consider the
// following tricky hp-scenario on two MPI processes:
//
//   subd 0 | subd 1
//          |
//  +---+---+---+---+
//  |   |   |   |   |
//  | 2 | 3 | 3 | 2 |
//  |   |   |   |   |
//  +---X-+-+-+-+---+
//  |   |2|2|2|2|   |
//  | 2 +-+-+-+-+ 2 |
//  |   |2|2|2|2|   |
//  +---+-+-+-+-+---+
//          |
//   subd 0 | subd 1
//
// Consider the edge marked with 'X'; in particular, DoFs 63 and 64 which solely
// belong to the Q3 element and are hanging nodes on that edge.
//
// The DoFs of the Q3 element are constrained against the coarse Q2 elements.
// The DoFs of the fine Q2 elements are also constrained against the coarse Q2
// elements. DoFTools::make_hanging_node_constraints() resolves these
// constraints correctly, and they are correct for process 0 on subdomain 0.
//
// From the point of view of process 1, the coarse Q2 elements on subdomain 0
// are artificial and their DoFs are thus not part of its locally relevant DoFs.
// Hence for the given context, DoFTools::make_hanging_node_constraints()
// constrains the DoFs of the Q3 element against the fine Q2 elements.
// As process 1 does not know about the coarse Q2 elements, i.e., for the given
// LOCAL context on process 1, these are indeed the correct constraints.
//
// However in the GLOBAL context, the constraints are wrong on process 1.
// They are only assembled correctly on process 0. Thus, we need
// AffineConstraints::make_consistent_in_parallel() to distribute these
// constraints correctly.


#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


void
test()
{
  constexpr int dim      = 3;
  constexpr int spacedim = 3;

  // set up grid
  parallel::distributed::Triangulation<dim, spacedim> triangulation(
    MPI_COMM_WORLD);
  {
    std::vector<unsigned int> repetitions({4, 2, 1});
    Point<dim>                bottom_left(-2, -1, 0);
    Point<dim>                top_right(2, 1, 1);

    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              repetitions,
                                              bottom_left,
                                              top_right);
  }

  // set up hp-refinement in center part
  DoFHandler<dim, spacedim> dof_handler(triangulation);
  {
    for (const auto &cell : dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      {
        const auto &center = cell->center();
        if (std::abs(center[0]) < 1.)
          {
            if (center[1] > 0.)
              cell->set_active_fe_index(1);
            else
              cell->set_refine_flag();
          }
      }
    triangulation.execute_coarsening_and_refinement();
  }

  // enumerate dofs
  hp::FECollection<dim, spacedim> fe_collection;
  {
    fe_collection.push_back(FE_Q<dim, spacedim>(2));
    fe_collection.push_back(FE_Q<dim, spacedim>(3));
    dof_handler.distribute_dofs(fe_collection);
  }

#if false
  // output vtu
  DataOut<dim, spacedim> data_out;
  {
    Vector<float> fe_degrees(triangulation.n_active_cells());
    for (const auto &cell : dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      fe_degrees(cell->active_cell_index()) = cell->get_fe().degree;

    Vector<float> subdomain(triangulation.n_active_cells());
    for (auto &subd : subdomain)
      subd = triangulation.locally_owned_subdomain();

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(fe_degrees, "fe_degrees");
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.build_patches();

    data_out.write_vtu_in_parallel("output.vtu",
                                   dof_handler.get_mpi_communicator());
  }
#endif

  // miscellaneous
  const auto show_constraints_63_64 =
    [&](const AffineConstraints<double> &constraints) {
      const unsigned int my_id =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      for (const unsigned int dof : {63, 64})
        {
          deallog << "What process " << my_id << " believes about DoF " << dof
                  << ":" << std::endl;
          for (const auto &c : constraints.get_lines())
            if (c.index == dof)
              for (const auto &entry : c.entries)
                deallog << "    constrained against " << entry.first
                        << " with weight " << entry.second << std::endl;
        }
    };

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet  locally_active_dofs =
    DoFTools::extract_locally_active_dofs(dof_handler);
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  // make hanging node constraints
  AffineConstraints<double> constraints;
  {
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

    deallog << "------------- make_hanging_node_constraints():" << std::endl;
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    show_constraints_63_64(constraints);
  }

  // make consistent in parallel on locally active dofs
  AffineConstraints<double> constraints_active;
  {
    constraints_active.copy_from(constraints);

    deallog << "------------- make_consistent_in_parallel() "
            << "on locally active dofs:" << std::endl;
    constraints_active.make_consistent_in_parallel(
      locally_owned_dofs,
      locally_active_dofs,
      dof_handler.get_mpi_communicator());
    show_constraints_63_64(constraints_active);
  }

  // make consistent in parallel on locally relevant dofs
  AffineConstraints<double> constraints_relevant;
  {
    constraints_relevant.copy_from(constraints);

    deallog << "------------- make_consistent_in_parallel() "
            << "on locally relevant dofs:" << std::endl;
    constraints_relevant.make_consistent_in_parallel(
      locally_owned_dofs,
      locally_relevant_dofs,
      dof_handler.get_mpi_communicator());
    show_constraints_63_64(constraints_relevant);
  }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  MPILogInitAll                    all;

  test();
}
