// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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


// Create a distributed triangulation with multigrid levels, copy it,
// and show that the CellIds have not changed.

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

using namespace dealii;

template <int dim>
void
test(int n_refinements, const int n_subdivisions, MPI_Comm comm)
{
  // create pdt
  parallel::distributed::Triangulation<dim> tria_pdt(
    comm,
    dealii::Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria_pdt, n_subdivisions);
  tria_pdt.refine_global(n_refinements);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form pdt
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria_pdt,
      comm,
      TriangulationDescription::Settings::construct_multigrid_hierarchy);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  // copare centers of cells
  for (auto &cell : tria_pft.active_cell_iterators())
    {
      CellId id        = cell->id();
      auto   cell_base = tria_pdt.create_cell_iterator(id);
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

  const int      dim  = 2;
  const MPI_Comm comm = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    const int n_refinements  = 3;
    const int n_subdivisions = 5;
    test<2>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }
  {
    deallog.push("3d");
    const int n_refinements  = 3;
    const int n_subdivisions = 4;
    test<3>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }
}
