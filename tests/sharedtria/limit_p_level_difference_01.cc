// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// verify restrictions on level differences imposed by
// hp::Refinement::limit_p_level_difference()
//
// set the center cell to the highest p-level in a hyper_cross geometry
// and verify that all other cells comply to the level difference


#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include <vector>

#include "../tests.h"


template <int dim>
void
test(const unsigned int fes_size,
     const unsigned int max_difference,
     const bool         allow_artificial_cells)
{
  Assert(fes_size > 0, ExcInternalError());
  Assert(max_difference > 0, ExcInternalError());

  // setup FE collection
  hp::FECollection<dim> fes;
  while (fes.size() < fes_size)
    fes.push_back(FE_Q<dim>(1));

  const unsigned int contains_fe_index = 0;
  const auto         sequence = fes.get_hierarchy_sequence(contains_fe_index);

  // setup cross-shaped mesh
  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD,
                                            Triangulation<dim>::none,
                                            allow_artificial_cells);
  {
    std::vector<unsigned int> sizes(Utilities::pow(2, dim),
                                    static_cast<unsigned int>(
                                      (sequence.size() - 1) / max_difference));
    GridGenerator::hyper_cross(tria, sizes);
  }

  deallog << "ncells:" << tria.n_cells() << ", nfes:" << fes.size() << std::endl
          << "sequence:" << sequence << std::endl;

  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fes);

  // set center cell to highest p-level
  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->is_locally_owned() && cell->center() == Point<dim>())
      cell->set_active_fe_index(sequence.back());

  const bool fe_indices_changed =
    hp::Refinement::limit_p_level_difference(dofh,
                                             max_difference,
                                             contains_fe_index);
  tria.execute_coarsening_and_refinement();

  (void)fe_indices_changed;
  Assert(fe_indices_changed, ExcInternalError());

  // display number of cells for each FE index
  std::vector<unsigned int> count(fes.size(), 0);
  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->is_locally_owned())
      count[cell->active_fe_index()]++;
  Utilities::MPI::sum(count, tria.get_communicator(), count);
  deallog << "fe count:" << count << std::endl;

#ifdef DEBUG
  // check each cell's active FE index by its distance from the center
  for (const auto &cell : dofh.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const double       distance = cell->center().distance(Point<dim>());
        const unsigned int expected_level =
          (sequence.size() - 1) -
          max_difference * static_cast<unsigned int>(std::round(distance));

        Assert(cell->active_fe_index() == sequence[expected_level],
               ExcInternalError());
      }
#endif

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
      test<2>(4, 1, allow_artificial_cells);
      test<2>(8, 2, allow_artificial_cells);
      test<2>(12, 3, allow_artificial_cells);
    }
}
