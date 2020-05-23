// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

/*
 * Purpose: checks compute_intergrid_transfer_representation
 * on distributed (p4est) triangulation
 *
 * Author: Alexander Grayver, 2015
 */

#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test(unsigned n_refinements)
{
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (rank == 0)
    {
      deallog << "Checking in " << dim << " space dimensions" << std::endl
              << "---------------------------------------" << std::endl;
    }

  // create serial grid and refine once
  parallel::distributed::Triangulation<dim> tria1(MPI_COMM_SELF);
  GridGenerator::hyper_cube(tria1, -1, 1);
  tria1.refine_global(1);

  // create distributed grid and refine twice
  parallel::distributed::Triangulation<dim> tria2(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria2, -1, 1);
  tria2.refine_global(n_refinements);

  // do some local refinement
  Point<dim> p0;
  p0 *= 0.;
  for (unsigned int i = 0; i < n_refinements; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator cell;
      for (cell = tria2.begin_active(); cell != tria2.end(); ++cell)
        {
          if (cell->is_locally_owned() &&
              (cell->center().distance(p0) < 0.71 / double(i + 1)))
            cell->set_refine_flag();
        }

      tria2.prepare_coarsening_and_refinement();
      tria2.execute_coarsening_and_refinement();
    }

  FE_DGQ<dim>     fe_dg(0);
  DoFHandler<dim> dof_handler1(tria1);
  dof_handler1.distribute_dofs(fe_dg);

  DoFHandler<dim> dof_handler2(tria2);
  dof_handler2.distribute_dofs(fe_dg);

  InterGridMap<DoFHandler<dim>> grid_1_to_2_map;
  grid_1_to_2_map.make_mapping(dof_handler1, dof_handler2);

  typedef std::vector<std::map<types::global_dof_index, float>> TransferRep;
  TransferRep transfer_representation;
  DoFTools::compute_intergrid_transfer_representation(
    dof_handler1, 0, dof_handler2, 0, grid_1_to_2_map, transfer_representation);

  // For this test case, all weights are one and their sum
  // should be equal to number of degrees of freedom
  unsigned local_sum = 0.;
  for (std::size_t i = 0; i < transfer_representation.size(); ++i)
    {
      TransferRep::value_type m = transfer_representation[i];
      for (TransferRep::value_type::const_iterator it = m.begin();
           it != m.end();
           ++it)
        local_sum += static_cast<unsigned int>(it->second);
    }

  unsigned global_sum = Utilities::MPI::sum(local_sum, MPI_COMM_WORLD);

  if (rank == 0)
    {
      deallog << "# dofs = " << dof_handler2.n_dofs() << std::endl;
      deallog << "sum(weights) = " << global_sum << std::endl;
      if (dof_handler2.n_dofs() == global_sum)
        deallog << "OK" << std::endl;
      else
        deallog << "Failed" << std::endl;
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(rank));

  unsigned n_refinements = 2;

  if (rank == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>(n_refinements);
      deallog.pop();
      deallog.push("3d");
      test<3>(n_refinements);
      deallog.pop();
    }
  else
    {
      test<2>(n_refinements);
      test<3>(n_refinements);
    }
}
