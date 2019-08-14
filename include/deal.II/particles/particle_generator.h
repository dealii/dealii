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

#ifndef dealii_particles_particle_generator_h
#define dealii_particles_particle_generator_h

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/particles/particle_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  /**
   * A namespace that contains all classes that are related to the particle
   * generation.
   */
  namespace Generator
  {
    /**
     * A function that generates particles in every cell at specified @p particle_reference_locations.
     */
    template <int dim, int spacedim = dim>
    void
    regular_reference_locations(
      const parallel::distributed::Triangulation<dim, spacedim> &triangulation,
      const std::vector<Point<dim>> & particle_reference_locations,
      ParticleHandler<dim, spacedim> &particle_handler,
      const Mapping<dim, spacedim> &  mapping = MappingQ1<dim, spacedim>());

  } // namespace Generator
} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif
