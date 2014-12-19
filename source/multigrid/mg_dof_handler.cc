// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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


#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/dofs/dof_faces.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/grid/tria_levels.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/exceptions.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    template <int, int> class Triangulation;
  }
}



/* ------------------------ MGDoFHandler ------------------------------------- */

template <int dim, int spacedim>
const unsigned int MGDoFHandler<dim,spacedim>::dimension;


template <int dim, int spacedim>
MGDoFHandler<dim,spacedim>::MGDoFHandler ()
{}



template <int dim, int spacedim>
MGDoFHandler<dim,spacedim>::MGDoFHandler (const Triangulation<dim,spacedim> &tria)
  :
  DoFHandler<dim,spacedim> (tria)
{}


template <int dim, int spacedim>
MGDoFHandler<dim,spacedim>::~MGDoFHandler ()
{}

template <int dim, int spacedim>
void MGDoFHandler<dim,spacedim>::distribute_dofs (const FiniteElement<dim,spacedim> &fe)
{
  // first distribute global dofs
  DoFHandler<dim,spacedim>::distribute_dofs (fe);
  DoFHandler<dim,spacedim>::distribute_mg_dofs (fe);
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::cell_iterator
MGDoFHandler<dim,spacedim>::begin (const unsigned int level) const
{
  return DoFHandler<dim,spacedim>::begin_mg(level);
}


template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::cell_iterator
MGDoFHandler<dim,spacedim>::end (const unsigned int level) const
{
  return DoFHandler<dim,spacedim>::end_mg(level);
}


template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::cell_iterator
MGDoFHandler<dim,spacedim>::end () const
{
  return DoFHandler<dim,spacedim>::end_mg();
}


// explicit instantiations
#include "mg_dof_handler.inst"


DEAL_II_NAMESPACE_CLOSE
