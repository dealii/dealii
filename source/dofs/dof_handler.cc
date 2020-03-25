// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
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

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler_base.templates.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
DoFHandler<dim, spacedim>::DoFHandler()
  : DoFHandlerBase<dim, spacedim, DoFHandler<dim, spacedim>>()
{}



template <int dim, int spacedim>
DoFHandler<dim, spacedim>::DoFHandler(const Triangulation<dim, spacedim> &tria)
  : DoFHandlerBase<dim, spacedim, DoFHandler<dim, spacedim>>(tria)
{}

/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_handler.inst"


DEAL_II_NAMESPACE_CLOSE
