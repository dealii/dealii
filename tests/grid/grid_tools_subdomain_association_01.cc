// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check GridTools::get_subdomain_association


#include <deal.II/base/mpi.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <memory>
#include <vector>

#include "../tests.h"

#include "../test_grids.h"


enum Type
{
  Shared,
  Distributed
};


template <int dim>
std::unique_ptr<parallel::TriangulationBase<dim>>
create_triangulation(Type type)
{
  if (type == Type::Shared)
    return std::make_unique<parallel::shared::Triangulation<dim>>(
      MPI_COMM_WORLD, Triangulation<dim>::none, true);
  else if (type == Type::Distributed)
    return std::make_unique<parallel::distributed::Triangulation<dim>>(
      MPI_COMM_WORLD);
  else
    return nullptr;
}


template <int dim>
void
test(parallel::TriangulationBase<dim> &tria)
{
  // ----- setup -----
  tria.clear();
  TestGrids::hyper_line(tria, 4);

  // refine the last cell
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->id().to_string() == "0_0:")
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  // ----- test -----
  // gather global cellid objects
  std::vector<CellId> global_cell_ids;
  global_cell_ids.reserve(tria.n_global_active_cells());
  {
    std::vector<CellId> local_cell_ids;
    local_cell_ids.reserve(tria.n_active_cells());
    for (const auto &cell :
         tria.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
      local_cell_ids.push_back(cell->id());

    std::vector<std::vector<CellId>> cell_ids_per_processor =
      Utilities::MPI::all_gather(MPI_COMM_WORLD, local_cell_ids);

    for (const auto &cell_ids : cell_ids_per_processor)
      global_cell_ids.insert(global_cell_ids.end(),
                             cell_ids.cbegin(),
                             cell_ids.cend());
  }

  // determine subdomain of every cellid
  std::vector<types::subdomain_id> subdomain_ids =
    GridTools::get_subdomain_association(tria, global_cell_ids);

  // ----- verify results -----
  AssertDimension(tria.n_global_active_cells(), global_cell_ids.size());
  AssertDimension(tria.n_global_active_cells(), subdomain_ids.size());

  // check if processor owns represented cells that have been assigned to them
  for (unsigned int i = 0; i < tria.n_global_active_cells(); ++i)
    if (tria.locally_owned_subdomain() == subdomain_ids[i])
      {
        const auto cell = tria.create_cell_iterator(global_cell_ids[i]);
        Assert(cell->is_locally_owned(), ExcInternalError());
      }

  // check if all processors have the same list of subdomain ids
  {
    std::vector<std::vector<types::subdomain_id>> subdomain_ids_per_processor =
      Utilities::MPI::gather(MPI_COMM_WORLD, subdomain_ids);

    if (tria.locally_owned_subdomain() == 0)
      for (unsigned int i = 1; i < subdomain_ids_per_processor.size(); ++i)
        Assert(subdomain_ids_per_processor[0] == subdomain_ids_per_processor[i],
               ExcInternalError());
  }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;


  for (const auto &type : {Type::Shared, Type::Distributed})
    {
      deallog.push("2d");
      test(*create_triangulation<2>(type));
      deallog.pop();
      deallog.push("3d");
      test(*create_triangulation<3>(type));
      deallog.pop();
    }
}
