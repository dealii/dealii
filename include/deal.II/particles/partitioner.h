// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_particles_partitioner_h
#define dealii_particles_partitioner_h

#include <deal.II/base/config.h>

#include <deal.II/particles/particle_iterator.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  namespace internal
  {
    /**
     * Cache structure used to store the elements which are required to
     * exchange the particle information (location and properties) across
     * processors in order to update the ghost particles.
     *
     * This structure should only be used when one wishes to carry out work
     * using the particles without calling
     * sort_particles_into_subdomain_and_cells at every iteration. This is
     * useful when particle-particle interaction occurs at a different time
     * scale than particle-mesh interaction.
     */
    template <int dim, int spacedim>
    struct GhostParticlePartitioner
    {
      /**
       * A type that can be used to iterate over all particles in the domain.
       */
      using particle_iterator = ParticleIterator<dim, spacedim>;

      /**
       * Indicates if the cache has been built to prevent updating particles
       * with an invalid cache.
       */
      bool valid = false;

      /**
       * Vector of the subdomain id of all possible neighbors of the current
       * subdomain.
       */
      std::vector<types::subdomain_id> neighbors;

      /**
       * Vector of size (neighbors.size()+1) used to store the start and the
       * end point of the data that must go from the current subdomain to the
       * neighbors. For neighbor i, send_pointers[i] indicates the beginning
       * and send_pointers[i+1] indicates the end of the data that must be
       * sent.
       */
      std::vector<unsigned int> send_pointers;

      /**
       * Set of particles that currently live in the ghost cells of the local
       * domain, organized by the subdomain_id. These
       * particles are equivalent to the ghost entries in distributed vectors.
       */
      std::map<types::subdomain_id, std::vector<particle_iterator>>
        ghost_particles_by_domain;

      /**
       * Vector of size (neighbors.size()+1) used to store the start and the
       * end point of the data that must be received from neighbor[i] on
       * the current subdomain. For neighbor i, recv_pointers[i] indicate the
       * beginning and recv_pointers[i+1] indicates the end of the data that
       * must be received.
       *
       * This structure is similar to
       * Utilities::MPI::Partitioner::import_targets when combined with
       * neighbors.
       */
      std::vector<unsigned int> recv_pointers;

      /**
       * Vector of ghost particles in the order in which they are inserted
       * in the multimap used to store particles on the triangulation. This
       * information is used to update the ghost particle information
       * without clearing the multimap of ghost particles, thus greatly
       * reducing the cost of exchanging the ghost particles information.
       */
      std::vector<particle_iterator> ghost_particles_iterators;

      /**
       * Temporary storage that holds the data of the particles to be sent
       * to other processors to update the ghost particles information
       * in update_ghost_particles()
       * send_recv_particles_properties_and_location()
       */
      std::vector<char> send_data;

      /**
       * Temporary storage that holds the data of the particles to receive
       * the ghost particles information from other processors in
       * update_ghost_particles()
       * send_recv_particles_properties_and_location()
       */
      std::vector<char> recv_data;
    };
  } // namespace internal

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif
