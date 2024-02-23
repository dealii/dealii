// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a serial triangulation with multigrid level, copy it, and
// show that the CellIds have not changed.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../grid/tests.h"


template <int dim>
void
test(int n_refinements, MPI_Comm comm)
{
  // create serial triangulation
  Triangulation<dim> basetria(
    Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_L(basetria);
  basetria.refine_global(n_refinements);

  GridTools::partition_triangulation_zorder(
    Utilities::MPI::n_mpi_processes(comm), basetria);
  GridTools::partition_multigrid_levels(basetria);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria,
      comm,
      TriangulationDescription::Settings::construct_multigrid_hierarchy);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);


  for (auto &cell : tria_pft.active_cell_iterators())
    {
      CellId id        = cell->id();
      auto   cell_base = basetria.create_cell_iterator(id);
      // Assert(cell->center() == cell_base->center(),
      //       ExcMessage("Cells do not match"));
      for (unsigned int d = 0; d < dim; ++d)
        Assert(std::abs(cell->center()[d] - cell_base->center()[d]) < 1e-9,
               ExcMessage("Cells do not match"));
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  const int      n_refinements = 3;
  const MPI_Comm comm          = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    test<2>(n_refinements, comm);
    deallog.pop();
  }
  {
    deallog.push("3d");
    test<3>(n_refinements, comm);
    deallog.pop();
  }
}
