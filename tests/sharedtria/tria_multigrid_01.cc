// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// create a shared tria mesh and partition
// active mesh in z-order as well as multigrid
// compare partition against p4est

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
void
compare_meshes(parallel::shared::Triangulation<dim>      &shared_tria,
               parallel::distributed::Triangulation<dim> &p4est_tria)
{
  std::map<CellId, unsigned int> shared_map;
  std::map<CellId, unsigned int> p4est_map;

  typename Triangulation<dim>::cell_iterator cell1 = shared_tria.begin(),
                                             endc1 = shared_tria.end();
  for (; cell1 != endc1; ++cell1)
    {
      if (cell1->is_locally_owned_on_level())
        shared_map.insert(
          std::make_pair(cell1->id(), cell1->level_subdomain_id()));
    }

  typename Triangulation<dim>::cell_iterator cell2 = p4est_tria.begin(),
                                             endc2 = p4est_tria.end();
  for (; cell2 != endc2; ++cell2)
    {
      if (cell2->is_locally_owned_on_level())
        p4est_map.insert(
          std::make_pair(cell2->id(), cell2->level_subdomain_id()));
    }

  for (std::map<CellId, unsigned int>::iterator it = p4est_map.begin();
       it != p4est_map.end();
       ++it)
    {
      AssertThrow(shared_map[it->first] == it->second,
                  ExcMessage(
                    "Not all CellIds map to correct level subdomain_ids."));
    }


  std::set<CellId> shared_known_cells;
  std::set<CellId> p4est_known_cells;

  cell1 = shared_tria.begin(), endc1 = shared_tria.end();
  for (; cell1 != endc1; ++cell1)
    {
      if (cell1->level_subdomain_id() != numbers::artificial_subdomain_id)
        shared_known_cells.insert(cell1->id());
    }
  cell2 = p4est_tria.begin(), endc2 = p4est_tria.end();
  for (; cell2 != endc2; ++cell2)
    {
      if (cell2->level_subdomain_id() != numbers::artificial_subdomain_id)
        p4est_known_cells.insert(cell2->id());
    }

  cell1 = shared_tria.begin(), endc1 = shared_tria.end();
  for (; cell1 != endc1; ++cell1)
    {
      if (cell1->level_subdomain_id() != numbers::artificial_subdomain_id)
        Assert(p4est_known_cells.find(cell1->id()) != p4est_known_cells.end(),
               ExcMessage("Cell not known by processor."));
    }
  cell2 = p4est_tria.begin(), endc2 = p4est_tria.end();
  for (; cell2 != endc2; ++cell2)
    {
      if (cell2->level_subdomain_id() != numbers::artificial_subdomain_id)
        Assert(shared_known_cells.find(cell2->id()) != shared_known_cells.end(),
               ExcMessage("Cell not known by processor."));
    }
}


template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> shared_tria(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::limit_level_difference_at_vertices),
    true,
    typename parallel::shared::Triangulation<dim>::Settings(
      parallel::shared::Triangulation<dim>::partition_zorder |
      parallel::shared::Triangulation<dim>::construct_multigrid_hierarchy));


  parallel::distributed::Triangulation<dim> p4est_tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);


  unsigned int refinements = 2;
  GridGenerator::subdivided_hyper_cube(shared_tria, 2, -1, 1);
  shared_tria.refine_global(refinements);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         shared_tria.begin_active();
       cell != shared_tria.end();
       ++cell)
    if (cell->center().norm() < 0.55)
      cell->set_refine_flag();
  shared_tria.execute_coarsening_and_refinement();
  for (typename Triangulation<dim>::active_cell_iterator cell =
         shared_tria.begin_active();
       cell != shared_tria.end();
       ++cell)
    if (cell->center().norm() > 0.3 && cell->center().norm() < 0.42)
      cell->set_refine_flag();
  shared_tria.execute_coarsening_and_refinement();
  for (typename Triangulation<dim>::active_cell_iterator cell =
         shared_tria.begin_active();
       cell != shared_tria.end();
       ++cell)
    if (cell->at_boundary() && (cell->center()[0] < 0 || cell->center()[1] < 0))
      cell->set_refine_flag();
  shared_tria.execute_coarsening_and_refinement();

  GridGenerator::subdivided_hyper_cube(p4est_tria, 2, -1, 1);
  p4est_tria.refine_global(refinements);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         p4est_tria.begin_active();
       cell != p4est_tria.end();
       ++cell)
    if (cell->is_locally_owned() && cell->center().norm() < 0.55)
      cell->set_refine_flag();
  p4est_tria.execute_coarsening_and_refinement();
  for (typename Triangulation<dim>::active_cell_iterator cell =
         p4est_tria.begin_active();
       cell != p4est_tria.end();
       ++cell)
    if (cell->is_locally_owned() && cell->center().norm() > 0.3 &&
        cell->center().norm() < 0.42)
      cell->set_refine_flag();
  p4est_tria.execute_coarsening_and_refinement();
  for (typename Triangulation<dim>::active_cell_iterator cell =
         p4est_tria.begin_active();
       cell != p4est_tria.end();
       ++cell)
    if (cell->is_locally_owned() && cell->at_boundary() &&
        (cell->center()[0] < 0 || cell->center()[1] < 0))
      cell->set_refine_flag();
  p4est_tria.execute_coarsening_and_refinement();

  compare_meshes(shared_tria, p4est_tria);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();

  deallog << "OK" << std::endl;
}
