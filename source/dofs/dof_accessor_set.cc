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
template <class OutputVector, typename number>
void
DoFCellAccessor<DH,lda>::
set_dof_values_by_interpolation (const Vector<number> &local_values,
                                 OutputVector         &values,
                                 const unsigned int fe_index) const
{
  if (!this->has_children())
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
        // simply set the values on this cell
        this->set_dof_values (local_values, values);
      else
        {
          Assert (local_values.size() == this->dof_handler->get_fe()[fe_index].dofs_per_cell,
                  ExcMessage ("Incorrect size of local_values vector.") );

          FullMatrix<double> interpolation (this->get_fe().dofs_per_cell, this->dof_handler->get_fe()[fe_index].dofs_per_cell);

          this->get_fe().get_interpolation_matrix (this->dof_handler->get_fe()[fe_index],
                                                   interpolation);

          // do the interpolation to the target space. for historical
          // reasons, matrices are set to size 0x0 internally even
          // we reinit as 4x0, so we have to treat this case specially
          Vector<number> tmp (this->get_fe().dofs_per_cell);
          if ((tmp.size() > 0) && (local_values.size() > 0))
            interpolation.vmult (tmp, local_values);

          // now set the dof values in the global vector
          this->set_dof_values (tmp, values);
        }
    }
  else
    // otherwise distribute them to the children
    {
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
      Assert (&this->get_dof_handler().get_fe() != 0,
              typename BaseClass::ExcInvalidObject());
      Assert (local_values.size() == dofs_per_cell,
              typename BaseClass::ExcVectorDoesNotMatch());
      Assert (values.size() == this->dof_handler->n_dofs(),
              typename BaseClass::ExcVectorDoesNotMatch());

      Vector<number> tmp(dofs_per_cell);

      for (unsigned int child=0; child<this->n_children(); ++child)
        {
          if (tmp.size() > 0)
            fe.get_prolongation_matrix(child, this->refinement_case())
            .vmult (tmp, local_values);
          this->child(child)->set_dof_values_by_interpolation (tmp, values, fe_index);
        }
    }
}


// --------------------------------------------------------------------------
// explicit instantiations
#include "dof_accessor_set.inst"

DEAL_II_NAMESPACE_CLOSE
