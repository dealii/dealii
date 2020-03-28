// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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

#include <deal.II/hp/dof_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  template <int dim, int spacedim>
  DoFHandler<dim, spacedim>::DoFHandler()
    : dealii::DoFHandler<dim, spacedim>(true)
  {}



  template <int dim, int spacedim>
  DoFHandler<dim, spacedim>::DoFHandler(
    const Triangulation<dim, spacedim> &tria)
    : dealii::DoFHandler<dim, spacedim>(tria, true)
  {}

} // namespace hp

/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_handler.inst"


DEAL_II_NAMESPACE_CLOSE
