// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// verify restrictions on level differences imposed by
// hp::Refinement::limit_p_level_difference()
// on h-refined grids


#include <deal.II/base/point.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include "../tests.h"

#include "../test_grids.h"


unsigned int
divide_and_ceil(unsigned int x, unsigned int y)
{
  return x / y + (x % y > 0);
}


template <int dim>
void
test(const unsigned int fes_size,
     const unsigned int max_difference,
     const bool         allow_artificial_cells)
{
  Assert(max_difference > 0, ExcInternalError());
  Assert(fes_size > 1, ExcInternalError());

  // setup FE collection
  hp::FECollection<dim> fes;
  while (fes.size() < fes_size)
    fes.push_back(FE_Q<dim>(1));

  const unsigned int contains_fe_index = 0;
  const auto         sequence = fes.get_hierarchy_sequence(contains_fe_index);

  // set up line grid
  // - refine leftmost column of cells consecutively
  // - assign highest p-level to rightmost cell
  //
  // +++-+---+------+
  // +++ |   |      |
  // +++-+---+      |
  // +++ |   |      |
  // +++-+---+------+
  //
  // after prepare_coarsening_and_refinement(), each p-level will correspond to
  // a unique column of cells and thus h-level

  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD,
                                            Triangulation<dim>::none,
                                            allow_artificial_cells);
  TestGrids::hyper_line(tria, 2);
  const unsigned int n_refinements =
    divide_and_ceil(sequence.size() - 1, max_difference);
  for (unsigned int i = 0; i < n_refinements; ++i)
    {
      for (const auto &cell : tria.active_cell_iterators())
        for (unsigned int v = 0; v < cell->n_vertices(); ++v)
          if (cell->vertex(v)[0] == 0.)
            cell->set_refine_flag();

      tria.execute_coarsening_and_refinement();
    }

  DoFHandler<dim> dofh(tria);
  for (const auto &cell : dofh.cell_iterators_on_level(0))
    if (cell->is_active() && cell->is_locally_owned())
      cell->set_active_fe_index(sequence.back());
  dofh.distribute_dofs(fes);

  const bool fe_indices_changed =
    hp::Refinement::limit_p_level_difference(dofh,
                                             max_difference,
                                             contains_fe_index);
  tria.execute_coarsening_and_refinement();

  (void)fe_indices_changed;
  Assert(fe_indices_changed, ExcInternalError());

  if constexpr (running_in_debug_mode())
    {
      // check each cell's active FE by its h-level
      for (unsigned int l = 0; l < tria.n_global_levels(); ++l)
        for (const auto &cell : dofh.cell_iterators_on_level(l))
          if (cell->is_active() && cell->is_locally_owned())
            {
              const unsigned int expected_level = std::max(
                0, static_cast<int>(sequence.size() - 1 - l * max_difference));
              Assert(cell->active_fe_index() == sequence[expected_level],
                     ExcInternalError());
            }
    }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    initlog();

  deallog << std::boolalpha;
  for (const bool allow_artificial_cells : {false, true})
    {
      deallog << "artificial cells: " << allow_artificial_cells << std::endl;
      test<2>(5, 1, allow_artificial_cells);
      test<2>(10, 2, allow_artificial_cells);
      test<2>(15, 3, allow_artificial_cells);
    }
}
