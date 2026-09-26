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


// Create a distributed triangulation without multigrid levels and copy it.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include "../grid/tests.h"

template <int dim>
void
distorte_coarse_grid(Triangulation<dim> &tria)
{
  const auto [points, cells, subcelldata] =
    GridTools::get_coarse_mesh_description(tria);

  std::vector<CellData<dim>> new_cells;

  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = i; j < cells.size(); j += 2)
      new_cells.emplace_back(cells[j]);

  tria.clear();
  tria.create_triangulation(points, new_cells, {});
}


template <int dim>
void
test(int n_refinements, MPI_Comm comm)
{
  // create pdt
  parallel::distributed::Triangulation<dim> tria_pdt(comm);
  GridGenerator::subdivided_hyper_cube(tria_pdt,
                                       Utilities::pow(2, n_refinements) + 1);
  distorte_coarse_grid(tria_pdt);

  deallog << "P::D::T:" << std::endl;
  for (const auto &cell : tria_pdt.active_cell_iterators())
    deallog << " " << cell->center() << " -> " << cell->subdomain_id()
            << std::endl;

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria_pdt, comm);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  deallog << "P::F::T:" << std::endl;
  for (const auto &cell : tria_pft.active_cell_iterators())
    deallog << " " << cell->center() << " -> " << cell->subdomain_id()
            << std::endl;


  // actually create triangulation
  construction_data.reorder_coarse_grid();
  tria_pft.clear();
  tria_pft.create_triangulation(construction_data);

  deallog << "P::F::T:" << std::endl;
  for (const auto &cell : tria_pft.active_cell_iterators())
    deallog << " " << cell->center() << " -> " << cell->subdomain_id()
            << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    const int n_refinements = 3;
    test<2>(n_refinements, comm);
    deallog.pop();
  }
}
