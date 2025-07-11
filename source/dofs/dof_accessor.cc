// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_levels.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria_iterator.h>

#include <boost/container/small_vector.hpp>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/* --------------------- Static variables: DoFAccessor --------------------- */

template <int structdim, int dim, int spacedim, bool level_dof_access>
const unsigned int
  DoFAccessor<structdim, dim, spacedim, level_dof_access>::dimension;

template <int structdim, int dim, int spacedim, bool level_dof_access>
const unsigned int
  DoFAccessor<structdim, dim, spacedim, level_dof_access>::space_dimension;



/* ---------------------- Functions: DoFInvalidAccessor -------------------- */

template <int structdim, int dim, int spacedim>
DoFInvalidAccessor<structdim, dim, spacedim>::DoFInvalidAccessor(
  const void *,
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
  const types::fe_index) const
{
  DEAL_II_ASSERT_UNREACHABLE();
  return 0;
}



template <int structdim, int dim, int spacedim>
void
DoFInvalidAccessor<structdim, dim, spacedim>::set_dof_index(
  const unsigned int,
  const types::global_dof_index,
  const types::fe_index) const
{
  DEAL_II_ASSERT_UNREACHABLE();
}



/*------------------------- Functions: DoFCellAccessor -----------------------*/



template <int dim, int spacedim, bool lda>
void
DoFCellAccessor<dim, spacedim, lda>::set_dof_indices(
  const std::vector<types::global_dof_index> &local_dof_indices)
{
  Assert(static_cast<unsigned int>(this->present_level) <
           this->dof_handler->object_dof_indices.size(),
         ExcMessage("DoFHandler not initialized"));

  Assert(this->dof_handler != nullptr, typename BaseClass::ExcInvalidObject());

  internal::DoFAccessorImplementation::Implementation::
    template set_dof_indices<dim, spacedim, lda, dim>(*this,
                                                      local_dof_indices,
                                                      this->active_fe_index());
}



template <int dim, int spacedim, bool lda>
TriaIterator<DoFCellAccessor<dim, spacedim, lda>>
DoFCellAccessor<dim, spacedim, lda>::neighbor_child_on_subface(
  const unsigned int face,
  const unsigned int subface) const
{
  const TriaIterator<CellAccessor<dim, spacedim>> q =
    CellAccessor<dim, spacedim>::neighbor_child_on_subface(face, subface);
  return TriaIterator<DoFCellAccessor<dim, spacedim, lda>>(*q,
                                                           this->dof_handler);
}



template <int dim, int spacedim, bool lda>
TriaIterator<DoFCellAccessor<dim, spacedim, lda>>
DoFCellAccessor<dim, spacedim, lda>::periodic_neighbor_child_on_subface(
  const unsigned int face,
  const unsigned int subface) const
{
  const TriaIterator<CellAccessor<dim, spacedim>> q =
    CellAccessor<dim, spacedim>::periodic_neighbor_child_on_subface(face,
                                                                    subface);
  return TriaIterator<DoFCellAccessor<dim, spacedim, lda>>(*q,
                                                           this->dof_handler);
}



template <int dim, int spacedim, bool lda>
inline TriaIterator<DoFCellAccessor<dim, spacedim, lda>>
DoFCellAccessor<dim, spacedim, lda>::periodic_neighbor(
  const unsigned int face) const
{
  const TriaIterator<CellAccessor<dim, spacedim>> q =
    CellAccessor<dim, spacedim>::periodic_neighbor(face);
  return TriaIterator<DoFCellAccessor<dim, spacedim, lda>>(*q,
                                                           this->dof_handler);
}



template <int dim, int spacedim, bool lda>
inline TriaIterator<DoFCellAccessor<dim, spacedim, lda>>
DoFCellAccessor<dim, spacedim, lda>::neighbor_or_periodic_neighbor(
  const unsigned int face) const
{
  const TriaIterator<CellAccessor<dim, spacedim>> q =
    CellAccessor<dim, spacedim>::neighbor_or_periodic_neighbor(face);
  return TriaIterator<DoFCellAccessor<dim, spacedim, lda>>(*q,
                                                           this->dof_handler);
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::get_dof_indices(
  std::vector<types::global_dof_index> &dof_indices,
  const types::fe_index                 fe_index_) const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());

  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  Assert(static_cast<unsigned int>(this->level()) <
           this->dof_handler->object_dof_indices.size(),
         ExcMessage(
           "The DoFHandler to which this accessor points has not "
           "been initialized, i.e., it doesn't appear that DoF indices "
           "have been distributed on it."));

  // this function really only makes sense if either a) there are degrees of
  // freedom defined on the present object, or b) the object is non-active
  // objects but all degrees of freedom are located on vertices, since
  // otherwise there are degrees of freedom on sub-objects which are not
  // allocated for this non-active thing
  Assert(this->fe_index_is_active(fe_index) ||
           (this->dof_handler->get_fe(fe_index).n_dofs_per_cell() ==
            this->n_vertices() *
              this->dof_handler->get_fe(fe_index).n_dofs_per_vertex()),
         ExcInternalError());

  // now do the actual work
  dealii::internal::DoFAccessorImplementation::Implementation::get_dof_indices(
    *this, dof_indices, fe_index);
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::get_mg_dof_indices(
  const int                             level,
  std::vector<types::global_dof_index> &dof_indices,
  const types::fe_index                 fe_index_) const
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());
  Assert(this->dof_handler->mg_vertex_dofs.size() > 0,
         ExcMessage("Multigrid DoF indices can only be accessed after "
                    "DoFHandler::distribute_mg_dofs() has been called!"));

  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  internal::DoFAccessorImplementation::Implementation::get_mg_dof_indices(
    *this, level, dof_indices, fe_index);
}



template <int structdim, int dim, int spacedim, bool level_dof_access>
inline void
DoFAccessor<structdim, dim, spacedim, level_dof_access>::set_mg_dof_indices(
  const int                                   level,
  const std::vector<types::global_dof_index> &dof_indices,
  const types::fe_index                       fe_index_)
{
  Assert(this->dof_handler != nullptr, ExcInvalidObject());

  const auto fe_index =
    internal::DoFAccessorImplementation::get_fe_index_or_default(*this,
                                                                 fe_index_);

  internal::DoFAccessorImplementation::Implementation::set_mg_dof_indices(
    *this, level, dof_indices, fe_index);
}



namespace internal
{
  namespace DoFAccessorImplementation
  {
    template <int dim, int spacedim, bool level_dof_access>
    void
    get_cell_dof_indices(
      const dealii::DoFCellAccessor<dim, spacedim, level_dof_access> &accessor,
      boost::container::small_vector<types::global_dof_index, 27> &dof_indices,
      const unsigned int                                           fe_index)
    {
      Implementation::process_dof_indices(
        accessor,
        dof_indices,
        fe_index,
        Implementation::DoFIndexProcessor<dim, spacedim>(),
        [](auto stored_index, auto dof_ptr) { *dof_ptr = stored_index; },
        false);
    }
  } // namespace DoFAccessorImplementation
} // namespace internal



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  get_dof_indices(std::vector<types::global_dof_index> &dof_indices) const
{
  Assert(this->is_active(),
         ExcMessage("get_dof_indices() only works on active cells."));
  Assert(this->is_artificial() == false,
         ExcMessage("Can't ask for DoF indices on artificial cells."));
  AssertDimension(dof_indices.size(), this->get_fe().n_dofs_per_cell());

  dealii::internal::DoFAccessorImplementation::Implementation::get_dof_indices(
    *this, dof_indices, this->active_fe_index());
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  get_mg_dof_indices(std::vector<types::global_dof_index> &dof_indices) const
{
  Assert(this->dof_handler->mg_vertex_dofs.size() > 0,
         ExcMessage("Multigrid DoF indices can only be accessed after "
                    "DoFHandler::distribute_mg_dofs() has been called!"));
  DoFAccessor<dimension_, dimension_, space_dimension_, level_dof_access>::
    get_mg_dof_indices(this->level(), dof_indices);
}



template <int dimension_, int space_dimension_, bool level_dof_access>
inline void
DoFCellAccessor<dimension_, space_dimension_, level_dof_access>::
  set_mg_dof_indices(const std::vector<types::global_dof_index> &dof_indices)
{
  Assert(this->dof_handler->mg_vertex_dofs.size() > 0,
         ExcMessage("Multigrid DoF indices can only be accessed after "
                    "DoFHandler::distribute_mg_dofs() has been called!"));
  DoFAccessor<dimension_, dimension_, space_dimension_, level_dof_access>::
    set_mg_dof_indices(this->level(), dof_indices);
}

// --------------------------------------------------------------------------
// explicit instantiations
#include "dofs/dof_accessor.inst"

DEAL_II_NAMESPACE_CLOSE
