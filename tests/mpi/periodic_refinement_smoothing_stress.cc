/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
 * Copyright (C) 2026 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * ------------------------------------------------------------------------
 */

// Stress periodic refinement and coarsening matching using
// deterministic random flags on larger 2d and 3d distributed meshes
// to test match_refinement_across_periodic_boundaries flag.
//
// You can run this program separately and specify the number of seconds to run
// as the first argument. By default it runs for 5 seconds.

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>

#include "../tests.h"


namespace
{
  std::uint64_t
  hash_string(const std::string &string)
  {
    std::uint64_t hash = 14695981039346656037ULL;
    for (const unsigned char character : string)
      {
        hash ^= character;
        hash *= 1099511628211ULL;
      }
    return hash;
  }



  std::uint64_t
  mix(std::uint64_t value)
  {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
  }



  template <int dim>
  void
  check_periodic_levels(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const unsigned int                               iteration)
  {
    unsigned int n_mismatches = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const unsigned int face : cell->face_indices())
          if (cell->has_periodic_neighbor(face) &&
              cell->level() != cell->periodic_neighbor(face)->level())
            {
              ++n_mismatches;
              std::cerr << "Rank "
                        << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                        << ", dim " << dim << ", iteration " << iteration
                        << ": periodic cell " << cell->id() << " at level "
                        << cell->level() << " has neighbor "
                        << cell->periodic_neighbor(face)->id() << " at level "
                        << cell->periodic_neighbor(face)->level() << std::endl;
            }

    AssertThrow(Utilities::MPI::sum(n_mismatches, MPI_COMM_WORLD) == 0,
                ExcMessage("Periodic neighbors have different levels."));
  }



  template <int dim>
  void
  run_for(parallel::distributed::Triangulation<dim>   &triangulation,
          const std::chrono::steady_clock::time_point &deadline,
          unsigned int                                &iteration)
  {
    bool continue_running;
    do
      {
        const bool allow_refinement =
          triangulation.n_global_active_cells() < 100000;
        const std::uint64_t iteration_seed =
          static_cast<std::uint64_t>(iteration) * 0x9e3779b97f4a7c15ULL;

        // First choose individual refinements. Coarsening is handled in a
        // second pass so that we can reject a coarsening request for an
        // entire sibling family if any of its children is being refined.
        for (const auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              const std::uint64_t cell_hash =
                hash_string(cell->id().to_string());
              const std::uint64_t refine_choice =
                mix(cell_hash ^ iteration_seed) % 100;

              if (allow_refinement && cell->level() < 7 && refine_choice < 5)
                cell->set_refine_flag();
            }

        for (const auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned() && cell->level() > 2)
            {
              const auto parent                  = cell->parent();
              bool       all_children_are_active = true;
              bool       any_child_will_refine   = false;
              for (unsigned int child = 0; child < parent->n_children();
                   ++child)
                {
                  const auto sibling = parent->child(child);
                  all_children_are_active &= sibling->is_active();
                  if (allow_refinement && sibling->level() < 7)
                    {
                      const std::uint64_t sibling_hash =
                        hash_string(sibling->id().to_string());
                      any_child_will_refine |=
                        mix(sibling_hash ^ iteration_seed) % 100 < 5;
                    }
                }

              const std::uint64_t parent_hash =
                hash_string(parent->id().to_string());
              const std::uint64_t coarsen_choice =
                mix(parent_hash ^ iteration_seed ^ 0xd1b54a32d192ed03ULL) % 100;
              if (all_children_are_active && !any_child_will_refine &&
                  coarsen_choice < 20)
                cell->set_coarsen_flag();
            }

        triangulation.execute_coarsening_and_refinement();
        std::cout << "cells: " << triangulation.n_global_active_cells()
                  << ", levels: " << triangulation.n_global_levels()
                  << std::endl;
        check_periodic_levels(triangulation, iteration);
        ++iteration;

        continue_running = Utilities::MPI::logical_or(
          std::chrono::steady_clock::now() < deadline, MPI_COMM_WORLD);
      }
    while (continue_running);
  }



  template <int dim>
  void
  test(const std::chrono::steady_clock::time_point &deadline)
  {
    using Smoothing = typename Triangulation<dim>::MeshSmoothing;
    parallel::distributed::Triangulation<dim> triangulation(
      MPI_COMM_WORLD,
      Smoothing(
        Triangulation<dim>::limit_level_difference_at_vertices |
        Triangulation<dim>::smoothing_on_refinement |
        Triangulation<dim>::smoothing_on_coarsening |
        Triangulation<dim>::match_refinement_across_periodic_boundaries));

    GridGenerator::hyper_cube(triangulation, 0., 1., true);

    std::vector<
      GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
      periodic_faces;
    GridTools::collect_periodic_faces(
      triangulation, 0, 1, /*direction=*/0, periodic_faces);
    triangulation.add_periodicity(periodic_faces);

    triangulation.refine_global(3);

    unsigned int current_iteration = 0;
    run_for(triangulation, deadline, current_iteration);
  }
} // namespace



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  mpi_initlog();

  // The optional argument sets the total run time in seconds. The default
  // keeps the regular test suite short; larger values are useful for stress
  // runs.
  const unsigned int run_time_in_seconds = argc > 1 ? std::stoi(argv[1]) : 5;
  const auto         start               = std::chrono::steady_clock::now();

  test<2>(start + std::chrono::seconds(run_time_in_seconds / 2));
  test<3>(start + std::chrono::seconds(run_time_in_seconds / 2));

  deallog << "OK" << std::endl;
}
