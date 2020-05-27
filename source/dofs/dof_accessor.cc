// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_levels.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_iterator.templates.h>

#include <deal.II/hp/dof_handler.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/* --------------------- Static variables: DoFAccessor --------------------- */

template <int structdim, typename DoFHandlerType, bool level_dof_access>
const unsigned int
  DoFAccessor<structdim, DoFHandlerType, level_dof_access>::dimension;

template <int structdim, typename DoFHandlerType, bool level_dof_access>
const unsigned int
  DoFAccessor<structdim, DoFHandlerType, level_dof_access>::space_dimension;



/* ---------------------- Functions: DoFInvalidAccessor -------------------- */

template <int structdim, int dim, int spacedim>
DoFInvalidAccessor<structdim, dim, spacedim>::DoFInvalidAccessor(
  const Triangulation<dim, spacedim> *,
  const int,
  const int,
  const AccessorData *)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
DoFInvalidAccessor<structdim, dim, spacedim>::DoFInvalidAccessor(
  const DoFInvalidAccessor &i)
  : InvalidAccessor<structdim, dim, spacedim>(
      static_cast<const InvalidAccessor<structdim, dim, spacedim> &>(i))
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
types::global_dof_index
DoFInvalidAccessor<structdim, dim, spacedim>::dof_index(
  const unsigned int,
  const unsigned int) const
{
  Assert(false, ExcInternalError());
  return 0;
}



template <int structdim, int dim, int spacedim>
void
DoFInvalidAccessor<structdim, dim, spacedim>::set_dof_index(
  const unsigned int,
  const types::global_dof_index,
  const unsigned int) const
{
  Assert(false, ExcInternalError());
}



/*------------------------- Functions: DoFCellAccessor -----------------------*/



template <typename DoFHandlerType, bool lda>
void
DoFCellAccessor<DoFHandlerType, lda>::update_cell_dof_indices_cache() const
{
  Assert(static_cast<unsigned int>(this->present_level) <
           this->dof_handler->levels.size(),
         ExcMessage("DoFHandler not initialized"));

  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());

  internal::DoFCellAccessorImplementation::Implementation::
    update_cell_dof_indices_cache(*this);
}



template <typename DoFHandlerType, bool lda>
void
DoFCellAccessor<DoFHandlerType, lda>::set_dof_indices(
  const std::vector<types::global_dof_index> &local_dof_indices)
{
  Assert(static_cast<unsigned int>(this->present_level) <
           this->dof_handler->levels.size(),
         ExcMessage("DoFHandler not initialized"));

  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());

  AssertDimension(local_dof_indices.size(), this->get_fe().dofs_per_cell);

  internal::DoFAccessorImplementation::Implementation::
    template set_dof_indices<DoFHandlerType, lda, DoFHandlerType::dimension>(
      *this, local_dof_indices, this->active_fe_index());
}



template <typename DoFHandlerType, bool lda>
TriaIterator<DoFCellAccessor<DoFHandlerType, lda>>
DoFCellAccessor<DoFHandlerType, lda>::neighbor_child_on_subface(
  const unsigned int face,
  const unsigned int subface) const
{
  const TriaIterator<CellAccessor<dim, spacedim>> q =
    CellAccessor<dim, spacedim>::neighbor_child_on_subface(face, subface);
  return TriaIterator<DoFCellAccessor<DoFHandlerType, lda>>(*q,
                                                            this->dof_handler);
}



template <typename DoFHandlerType, bool lda>
TriaIterator<DoFCellAccessor<DoFHandlerType, lda>>
DoFCellAccessor<DoFHandlerType, lda>::periodic_neighbor_child_on_subface(
  const unsigned int face,
  const unsigned int subface) const
{
  const TriaIterator<CellAccessor<dim, spacedim>> q =
    CellAccessor<dim, spacedim>::periodic_neighbor_child_on_subface(face,
                                                                    subface);
  return TriaIterator<DoFCellAccessor<DoFHandlerType, lda>>(*q,
                                                            this->dof_handler);
}



template <typename DoFHandlerType, bool lda>
inline TriaIterator<DoFCellAccessor<DoFHandlerType, lda>>
DoFCellAccessor<DoFHandlerType, lda>::periodic_neighbor(
  const unsigned int face) const
{
  const TriaIterator<CellAccessor<dim, spacedim>> q =
    CellAccessor<dim, spacedim>::periodic_neighbor(face);
  return TriaIterator<DoFCellAccessor<DoFHandlerType, lda>>(*q,
                                                            this->dof_handler);
}



template <typename DoFHandlerType, bool lda>
inline TriaIterator<DoFCellAccessor<DoFHandlerType, lda>>
DoFCellAccessor<DoFHandlerType, lda>::neighbor_or_periodic_neighbor(
  const unsigned int face) const
{
  const TriaIterator<CellAccessor<dim, spacedim>> q =
    CellAccessor<dim, spacedim>::neighbor_or_periodic_neighbor(face);
  return TriaIterator<DoFCellAccessor<DoFHandlerType, lda>>(*q,
                                                            this->dof_handler);
}



// --------------------------------------------------------------------------
// explicit instantiations
#include "dof_accessor.inst"

DEAL_II_NAMESPACE_CLOSE
