// ---------------------------------------------------------------------
// $Id$
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

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_iterator.templates.h>
#include <deal.II/fe/fe.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/*------------------------- Static variables: DoFAccessor -----------------------*/

template <int structdim, class DH, bool level_dof_access>
const unsigned int DoFAccessor<structdim,DH,level_dof_access>::dimension;

template <int structdim, class DH, bool level_dof_access>
const unsigned int DoFAccessor<structdim,DH,level_dof_access>::space_dimension;



/*------------------------- Functions: DoFCellAccessor -----------------------*/




template <class DH, bool lda>
void
DoFCellAccessor<DH,lda>::update_cell_dof_indices_cache () const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0, typename BaseClass::ExcInvalidObject());

  internal::DoFCellAccessor::Implementation::
  update_cell_dof_indices_cache (*this);
}



template <class DH, bool lda>
void
DoFCellAccessor<DH,lda>::set_dof_indices (const std::vector<types::global_dof_index> &local_dof_indices)
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0, typename BaseClass::ExcInvalidObject());

  internal::DoFCellAccessor::Implementation::
  set_dof_indices (*this, local_dof_indices);
}




template <class DH, bool lda>
TriaIterator<DoFCellAccessor<DH,lda> >
DoFCellAccessor<DH,lda>::neighbor_child_on_subface (const unsigned int face,
                                                    const unsigned int subface) const
{
  const TriaIterator<CellAccessor<dim,spacedim> > q
    = CellAccessor<dim,spacedim>::neighbor_child_on_subface (face, subface);
  return TriaIterator<DoFCellAccessor<DH,lda> > (*q, this->dof_handler);
}



template <class DH, bool lda>
template <class InputVector, typename number>
void
DoFCellAccessor<DH,lda>::
get_interpolated_dof_values (const InputVector &values,
                             Vector<number>    &interpolated_values) const
{
  const FiniteElement<dim,spacedim> &fe            = this->get_fe();
  const unsigned int        dofs_per_cell = fe.dofs_per_cell;

  Assert (this->dof_handler != 0,
          typename BaseClass::ExcInvalidObject());
  Assert (&fe != 0,
          typename BaseClass::ExcInvalidObject());
  Assert (interpolated_values.size() == dofs_per_cell,
          typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
          typename BaseClass::ExcVectorDoesNotMatch());

  if (!this->has_children())
    // if this cell has no children: simply return the exact values on this
    // cell
    this->get_dof_values (values, interpolated_values);
  else
    // otherwise clobber them from the children
    {
      Vector<number> tmp1(dofs_per_cell);
      Vector<number> tmp2(dofs_per_cell);

      interpolated_values = 0;

      // later on we will have to push the values interpolated from the child
      // to the mother cell into the output vector. unfortunately, there are
      // two types of elements: ones where you add up the contributions from
      // the different child cells, and ones where you overwrite.
      //
      // an example for the first is piecewise constant (and discontinuous)
      // elements, where we build the value on the coarse cell by averaging
      // the values from the cell (i.e. by adding up a fraction of the values
      // of their values)
      //
      // an example for the latter are the usual continuous elements. the
      // value on a vertex of a coarse cell must there be the same,
      // irrespective of the adjacent cell we are presently on. so we always
      // overwrite. in fact, we must, since we cannot know in advance how many
      // neighbors there will be, so there is no way to compute the average
      // with fixed factors
      //
      // so we have to find out to which type this element belongs. the
      // difficulty is: the finite element may be a composed one, so we can
      // only hope to do this for each shape function individually. in fact,
      // there are even weird finite elements (for example the Raviart-Thomas
      // element) which have shape functions that are additive (interior ones)
      // and others that are overwriting (face degrees of freedom that need to
      // be continuous across the face). to avoid checking this over and over
      // again, we do this once now and cache the results
      std::vector<bool> restriction_is_additive (dofs_per_cell);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        restriction_is_additive[i] = fe.restriction_is_additive(i);

      for (unsigned int child=0; child<this->n_children(); ++child)
        {
          // get the values from the present child, if necessary by
          // interpolation itself
          this->child(child)->get_interpolated_dof_values (values,
                                                           tmp1);
          // interpolate these to the mother cell
          fe.get_restriction_matrix(child, this->refinement_case()).vmult (tmp2, tmp1);

          // and add up or set them in the output vector
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            if (restriction_is_additive[i])
              interpolated_values(i) += tmp2(i);
            else if (tmp2(i) != number())
              interpolated_values(i) = tmp2(i);
        }
    }
}



template <class DH, bool lda>
template <class OutputVector, typename number>
void
DoFCellAccessor<DH,lda>::
set_dof_values_by_interpolation (const Vector<number> &local_values,
                                 OutputVector         &values) const
{
  const unsigned int dofs_per_cell = this->get_fe().dofs_per_cell;

  Assert (this->dof_handler != 0,
          typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0,
          typename BaseClass::ExcInvalidObject());
  Assert (local_values.size() == dofs_per_cell,
          typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
          typename BaseClass::ExcVectorDoesNotMatch());

  if (!this->has_children())
    // if this cell has no children: simply set the values on this cell
    this->set_dof_values (local_values, values);
  else
    // otherwise distribute them to the children
    {
      Vector<number> tmp(dofs_per_cell);

      for (unsigned int child=0; child<this->n_children(); ++child)
        {
          Assert (this->child(child)->get_fe().dofs_per_cell == dofs_per_cell,
                  ExcNotImplemented());

          // prolong the given data to the present cell. FullMatrix only wants
          // us to call vmult if the matrix size is actually non-zero, so
          // check that case
          if (tmp.size() > 0)
            {
              this->get_fe().get_prolongation_matrix(child, this->refinement_case())
              .vmult (tmp, local_values);

              this->child(child)->set_dof_values_by_interpolation (tmp, values);
            }
        }
    }
}



// --------------------------------------------------------------------------
// explicit instantiations
#include "dof_accessor.inst"

DEAL_II_NAMESPACE_CLOSE
