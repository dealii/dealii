// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/particle_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <int dim, int spacedim>
  void
  DataOut<dim, spacedim>::build_patches(
    const Particles::ParticleHandler<dim, spacedim> &particles)
  {
    dataset_names.clear();
    dataset_names.emplace_back("id");
    patches.resize(particles.n_locally_owned_particles());

    auto particle = particles.begin();
    for (unsigned int i = 0; particle != particles.end(); ++particle, ++i)
      {
        patches[i].vertices[0]    = particle->get_location();
        patches[i].patch_index    = i;
        patches[i].n_subdivisions = 1;
        patches[i].data.reinit(dataset_names.size(), 1);

        patches[i].data(0, 0) = particle->get_id();
      }
  }



  template <int dim, int spacedim>
  const std::vector<DataOutBase::Patch<0, spacedim>> &
  DataOut<dim, spacedim>::get_patches() const
  {
    return patches;
  }



  template <int dim, int spacedim>
  std::vector<std::string>
  DataOut<dim, spacedim>::get_dataset_names() const
  {
    return dataset_names;
  }

} // namespace Particles

#include "data_out.inst"

DEAL_II_NAMESPACE_CLOSE
