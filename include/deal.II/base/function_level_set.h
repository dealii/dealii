// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_function_level_set_h
#define dealii_function_level_set_h

#include <deal.II/base/config.h>

#include <deal.II/base/function_signed_distance.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  namespace LevelSet
  {
    template <int dim>
    using Sphere DEAL_II_DEPRECATED = SignedDistance::Sphere<dim>;

    template <int dim>
    using Plane DEAL_II_DEPRECATED = SignedDistance::Plane<dim>;

  } // namespace LevelSet
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
