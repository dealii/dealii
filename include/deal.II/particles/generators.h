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

#ifdef DEAL_II_WITH_P4EST

namespace Particles
{
  /**
   * A namespace that contains all classes that are related to the particle
   * generation.
   */
  namespace Generators
  {
    /**
     * A function that generates particles in every cell at specified @p particle_reference_locations.
     * The total number of particles that is added to the @p particle_handler object is
     * the number of locally owned cells of the @p triangulation times the number of
     * locations in @p particle_reference_locations. An optional @p mapping argument
     * can be used to map from @p particle_reference_locations to the real particle locations.
     *
     * @param triangulation The triangulation associated with the @p particle_handler.
     *
     * @param particle_reference_locations A vector of positions in the unit cell.
     * Particles will be generated in every cell at these locations.
     *
     * @param particle_handler The particle handler that will take ownership
     * of the generated particles.
     *
     * @param mapping An optional mapping object that is used to map reference
     * location in the unit cell to the real cells of the triangulation. If no
     * mapping is provided a MappingQ1 is assumed.
     */
    template <int dim, int spacedim = dim>
    void
    regular_reference_locations(
      const Triangulation<dim, spacedim> &triangulation,
      const std::vector<Point<dim>> &     particle_reference_locations,
      ParticleHandler<dim, spacedim> &    particle_handler,
      const Mapping<dim, spacedim> &      mapping =
        StaticMappingQ1<dim, spacedim>::mapping);

  } // namespace Generators
} // namespace Particles

#endif // DEAL_II_WITH_P4EST

DEAL_II_NAMESPACE_CLOSE

#endif
