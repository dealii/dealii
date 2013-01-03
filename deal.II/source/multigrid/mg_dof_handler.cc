//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

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
