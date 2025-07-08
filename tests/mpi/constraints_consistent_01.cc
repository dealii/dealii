// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test AffineConstraints<double>::is_consistent_in_parallel

#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
check(parallel::distributed::Triangulation<dim> &tria)
{
  MPILogInitAll all;

  DoFHandler<2> dof_handler(tria);
  FESystem<dim> fe(FE_Q<dim>(2), dim);

  dof_handler.distribute_dofs(fe);

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet  locally_active_dofs =
    DoFTools::extract_locally_active_dofs(dof_handler);
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  AffineConstraints<double> constraints;

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  for (unsigned int id = 0; id < 1; ++id)
    dealii::VectorTools::interpolate_boundary_values(
      dof_handler,
      id,
      Functions::ConstantFunction<dim>(static_cast<double>(id) + 42.0, dim),
      constraints);

  VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                  0, /*first component*/
                                                  std::set<types::boundary_id>{
                                                    2},
                                                  constraints);

  constraints.close();
  deallog << "LocallyOwned = " << std::flush;
  locally_owned_dofs.print(deallog.get_file_stream());
  deallog << "constraints:" << std::endl;
  constraints.print(deallog.get_file_stream());
  deallog << "consistent? "
          << constraints.is_consistent_in_parallel(
               Utilities::MPI::all_gather(MPI_COMM_WORLD,
                                          dof_handler.locally_owned_dofs()),
               locally_active_dofs,
               MPI_COMM_WORLD,
               true)
          << std::endl;
}



template <int dim>
void
test()
{
  {
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(tr, 0.0, 1.0, false);
    tr.refine_global(2);
    check(tr);
  }
  {
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(tr, 0.0, 1.0, true);
    tr.refine_global(2);
    check(tr);
  }
  {
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

    GridGenerator::subdivided_hyper_cube(tr, 2);
    tr.refine_global(1);
    for (auto cell = tr.begin_active(); cell != tr.end(); ++cell)
      {
        if (cell->id().to_string() == "0_1:0")
          cell->set_refine_flag();
        else if (cell->parent()->id().to_string() == "3_0:")
          cell->set_coarsen_flag();
      }

    tr.execute_coarsening_and_refinement();
    check(tr);
  }
  {
    deallog << "this should be inconsistent with more than one rank:"
            << std::endl;
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(tr, 0.0, 1.0, false);
    tr.refine_global(2);

    for (auto cell = tr.begin_active(); cell != tr.end(); ++cell)
      if (cell->is_locally_owned())
        {
          for (unsigned int f(0); f < GeometryInfo<2>::faces_per_cell; ++f)
            {
              if (cell->face(f)->at_boundary() &&
                  cell->face(f)->center()[0] < 1e-10)
                cell->face(f)->set_all_boundary_ids(1);
            }
        }

    check(tr);
  }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  test<2>();
}
