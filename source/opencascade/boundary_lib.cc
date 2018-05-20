// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/opencascade/boundary_lib.h>

#ifdef DEAL_II_WITH_OPENCASCADE

DEAL_II_NAMESPACE_OPEN

namespace OpenCASCADE
{
  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  NormalProjectionBoundary<dim, spacedim>::clone() const
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(
      new NormalProjectionBoundary(this->sh, this->tolerance));
  }

  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  DirectionalProjectionBoundary<dim, spacedim>::clone() const
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(
      new DirectionalProjectionBoundary(
        this->sh, this->direction, this->tolerance));
  }

  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  NormalToMeshProjectionBoundary<dim, spacedim>::clone() const
  {
    return std::unique_ptr<Manifold<dim, spacedim>>(
      new NormalToMeshProjectionBoundary(this->sh, this->tolerance));
  }
// Explicit instantiations
#  include "boundary_lib.inst"
} // namespace OpenCASCADE

DEAL_II_NAMESPACE_CLOSE

#endif
