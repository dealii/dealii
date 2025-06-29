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

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim, bool lda>
template <typename Number>
void
DoFCellAccessor<dim, spacedim, lda>::get_interpolated_dof_values(
  const ReadVector<Number> &values,
  Vector<Number>           &interpolated_values,
  const types::fe_index     fe_index_) const
{
  const types::fe_index fe_index =
    (this->dof_handler->hp_capability_enabled == false &&
     fe_index_ == numbers::invalid_fe_index) ?
      DoFHandler<dim, spacedim>::default_fe_index :
      fe_index_;

  if (this->is_active())
    // If this cell is active: simply return the exact values on this
    // cell unless the finite element we need to interpolate to is different
    // than the one we have on the current cell
    {
      if ((this->dof_handler->hp_capability_enabled == false) ||
          // for hp-DoFHandlers, we need to require that on
          // active cells, you either don't specify an fe_index,
          // or that you specify the correct one
          (fe_index == this->active_fe_index()) ||
          (fe_index == numbers::invalid_fe_index))
        this->get_dof_values(values, interpolated_values);
      else
        {
          // well, here we need to first get the values from the current
          // cell and then interpolate it to the element requested. this
          // can clearly only happen for DoFHandler objects in hp-mode
          const unsigned int dofs_per_cell = this->get_fe().n_dofs_per_cell();
          if (dofs_per_cell == 0)
            {
              interpolated_values = 0;
            }
          else
            {
              Vector<Number> tmp(dofs_per_cell);
              this->get_dof_values(values, tmp);

              FullMatrix<double> interpolation(
                this->dof_handler->get_fe(fe_index).n_dofs_per_cell(),
                this->get_fe().n_dofs_per_cell());
              this->dof_handler->get_fe(fe_index).get_interpolation_matrix(
                this->get_fe(), interpolation);
              interpolation.vmult(interpolated_values, tmp);
            }
        }
    }
  else
    // The cell is not active; we need to obtain data them from
    // children recursively.
    {
      // we are on a non-active cell. these do not have any finite
      // element associated with them in the hp-context (in the non-hp-
      // context, we can simply assume that the FE space to which we
      // want to interpolate is the same as for all elements in the
      // mesh). consequently, we cannot interpolate from children's FE
      // space to this cell's (unknown) FE space unless an explicit
      // fe_index is given
      Assert((this->dof_handler->hp_capability_enabled == false) ||
               (fe_index != numbers::invalid_fe_index),
             ExcMessage(
               "You cannot call this function on non-active cells "
               "of DoFHandler objects unless you provide an explicit "
               "finite element index because they do not have naturally "
               "associated finite element spaces associated: degrees "
               "of freedom are only distributed on active cells for which "
               "the active FE index has been set."));

      const FiniteElement<dim, spacedim> &fe =
        this->get_dof_handler().get_fe(fe_index);
      const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

      Assert(this->dof_handler != nullptr,
             typename BaseClass::ExcInvalidObject());
      Assert(interpolated_values.size() == dofs_per_cell,
             typename BaseClass::ExcVectorDoesNotMatch());
      Assert(values.size() == this->dof_handler->n_dofs(),
             typename BaseClass::ExcVectorDoesNotMatch());


      // see if the finite element we have on the current cell has any
      // degrees of freedom to begin with; if not (e.g., when
      // interpolating FE_Nothing), then simply skip all of the
      // following since the output vector would be of size zero
      // anyway (and in fact is of size zero, see the assertion above)
      if (fe.n_dofs_per_cell() > 0)
        {
          Vector<Number> tmp1(dofs_per_cell);
          Vector<Number> tmp2(dofs_per_cell);

          interpolated_values = 0;

          // later on we will have to push the values interpolated from the
          // child to the mother cell into the output vector. unfortunately,
          // there are two types of elements: ones where you add up the
          // contributions from the different child cells, and ones where you
          // overwrite.
          //
          // an example for the first is piecewise constant (and discontinuous)
          // elements, where we build the value on the coarse cell by averaging
          // the values from the cell (i.e. by adding up a fraction of the
          // values of their values)
          //
          // an example for the latter are the usual continuous elements. the
          // value on a vertex of a coarse cell must there be the same,
          // irrespective of the adjacent cell we are presently on. so we always
          // overwrite. in fact, we must, since we cannot know in advance how
          // many neighbors there will be, so there is no way to compute the
          // average with fixed factors
          //
          // so we have to find out to which type this element belongs. the
          // difficulty is: the finite element may be a composed one, so we can
          // only hope to do this for each shape function individually. in fact,
          // there are even weird finite elements (for example the
          // Raviart-Thomas element) which have shape functions that are
          // additive (interior ones) and others that are overwriting (face
          // degrees of freedom that need to be continuous across the face).
          for (unsigned int child = 0; child < this->n_children(); ++child)
            {
              // get the values from the present child, if necessary by
              // interpolation itself either from its own children or
              // by interpolating from the finite element on an active
              // child to the finite element space requested here
              this->child(child)->get_interpolated_dof_values(values,
                                                              tmp1,
                                                              fe_index);
              // interpolate these to the mother cell
              fe.get_restriction_matrix(child, this->refinement_case())
                .vmult(tmp2, tmp1);

              // and add up or set them in the output vector
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                if (fe.restriction_is_additive(i))
                  interpolated_values(i) += tmp2(i);
                else if (tmp2(i) != Number())
                  interpolated_values(i) = tmp2(i);
            }
        }
    }
}


// --------------------------------------------------------------------------
// explicit instantiations
#include "dofs/dof_accessor_get.inst"

DEAL_II_NAMESPACE_CLOSE
