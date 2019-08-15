// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#include <deal.II/base/signaling_nan.h>

#include <deal.II/particles/particle_generator.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  namespace Generator
  {
    template <int dim, int spacedim = dim>
    void
    regular_reference_locations(
      const parallel::distributed::Triangulation<dim, spacedim> &triangulation,
      const std::vector<Point<dim>> & particle_reference_locations,
      ParticleHandler<dim, spacedim> &particle_handler,
      const Mapping<dim, spacedim> &  mapping)
    {
      types::particle_index n_particles_to_generate =
        triangulation.n_locally_owned_active_cells() *
        particle_reference_locations.size();

      // The local start index is the number of all particles
      // generated on lower MPI ranks.
      types::particle_index local_start_index = 0;
      MPI_Scan(&n_particles_to_generate,
               &local_start_index,
               1,
               DEAL_II_PARTICLE_INDEX_MPI_TYPE,
               MPI_SUM,
               triangulation.get_communicator());
      local_start_index -= n_particles_to_generate;

      types::particle_index particle_index = local_start_index;

      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (const auto &reference_location :
                   particle_reference_locations)
                {
                  const Point<spacedim> position_real =
                    mapping.transform_unit_to_real_cell(cell,
                                                        reference_location);

                  const Particle<dim, spacedim> particle(position_real,
                                                         reference_location,
                                                         particle_index);
                  particle_handler.insert_particle(particle, cell);
                  ++particle_index;
                }
            }
        }

      particle_handler.update_cached_numbers();
    }
  } // namespace Generator
} // namespace Particles



DEAL_II_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_OPEN

#include "particle_generator.inst"

DEAL_II_NAMESPACE_CLOSE
