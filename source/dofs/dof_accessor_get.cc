// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
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


template <class DH, bool lda>
template <class InputVector, typename number>
void
DoFCellAccessor<DH,lda>::
get_interpolated_dof_values (const InputVector &values,
                             Vector<number>    &interpolated_values,
                             const unsigned int fe_index) const
{
  if (!this->has_children())
    // if this cell has no children: simply return the exact values on this
    // cell unless the finite element we need to interpolate to is different than
    // the one we have on the current cell
    {
      if ((dynamic_cast<DoFHandler<DH::dimension,DH::space_dimension>*>
           (this->dof_handler)
           != 0)
          ||
          // for hp-DoFHandlers, we need to require that on
          // active cells, you either don't specify an fe_index,
          // or that you specify the correct one
          (fe_index == this->active_fe_index())
          ||
          (fe_index == DH::default_fe_index))
        this->get_dof_values (values, interpolated_values);
      else
        {
          // well, here we need to first get the values from the current
          // cell and then interpolate it to the element requested. this
          // can clearly only happen for hp::DoFHandler objects
          Vector<number> tmp (this->get_fe().dofs_per_cell);
          this->get_dof_values (values, tmp);

          FullMatrix<double> interpolation (this->dof_handler->get_fe()[fe_index].dofs_per_cell,
                                            this->get_fe().dofs_per_cell);
          this->dof_handler->get_fe()[fe_index].get_interpolation_matrix (this->get_fe(),
              interpolation);
          interpolation.vmult (interpolated_values, tmp);
        }
    }
  else
    // otherwise obtain them from the children
    {
      // we are on a non-active cell. these do not have any finite
      // element associated with them in the hp context (in the non-hp
      // context, we can simply assume that the FE space to which we
      // want to interpolate is the same as for all elements in the
      // mesh). consequently, we cannot interpolate from children's FE
      // space to this cell's (unknown) FE space unless an explicit
      // fe_index is given
      Assert ((dynamic_cast<DoFHandler<DH::dimension,DH::space_dimension>*>
               (this->dof_handler)
               != 0)
              ||
              (fe_index != DH::default_fe_index),
              ExcMessage ("You cannot call this function on non-active cells "
                          "of hp::DoFHandler objects unless you provide an explicit "
                          "finite element index because they do not have naturally "
                          "associated finite element spaces associated: degrees "
                          "of freedom are only distributed on active cells for which "
                          "the active_fe_index has been set."));

      const FiniteElement<dim,spacedim> &fe            = this->get_dof_handler().get_fe()[fe_index];
      const unsigned int                 dofs_per_cell = fe.dofs_per_cell;

      Assert (this->dof_handler != 0,
              typename BaseClass::ExcInvalidObject());
      Assert (&fe != 0,
              typename BaseClass::ExcInvalidObject());
      Assert (interpolated_values.size() == dofs_per_cell,
              typename BaseClass::ExcVectorDoesNotMatch());
      Assert (values.size() == this->dof_handler->n_dofs(),
              typename BaseClass::ExcVectorDoesNotMatch());


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
          // interpolation itself either from its own children or
          // by interpolating from the finite element on an active
          // child to the finite element space requested here
          this->child(child)->get_interpolated_dof_values (values,
                                                           tmp1,
                                                           fe_index);
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


// --------------------------------------------------------------------------
// explicit instantiations
#include "dof_accessor_get.inst"

DEAL_II_NAMESPACE_CLOSE
