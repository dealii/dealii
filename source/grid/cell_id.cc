// ---------------------------------------------------------------------
//
// Copyright (C) 2015,2016 by the deal.II authors
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

#include <deal.II/grid/cell_id.h>

#include <deal.II/grid/tria.h>

#include <sstream>

DEAL_II_NAMESPACE_OPEN

std::string
CellId::to_string() const
{
  std::ostringstream ss;
  ss << *this;
  return ss.str();
}

template<int dim, int spacedim>
typename Triangulation<dim,spacedim>::cell_iterator
CellId::to_cell(const Triangulation<dim,spacedim> &tria) const
{
  typename Triangulation<dim,spacedim>::cell_iterator cell (&tria,0,coarse_cell_id);

  for (unsigned int i = 0; i < child_indices.size(); ++i)
    cell = cell->child(static_cast<unsigned int> (child_indices[i]));

  return cell;
}

// explicit instantiations
#include "cell_id.inst"

DEAL_II_NAMESPACE_CLOSE
