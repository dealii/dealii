//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>
#include <lac/sparse_matrix.h>

#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_levels.h>
#include <hp/dof_handler.h>
#include <grid/tria_boundary.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/*------------------------- Functions: DoFCellAccessor -----------------------*/




template <class DH>
void
DoFCellAccessor<DH>::update_cell_dof_indices_cache () const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0, typename BaseClass::ExcInvalidObject());

  internal::DoFCellAccessor::Implementation::
    update_cell_dof_indices_cache (*this);
}



template <class DH>
void
DoFCellAccessor<DH>::set_dof_indices (const std::vector<unsigned int> &local_dof_indices)
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0, typename BaseClass::ExcInvalidObject());

  internal::DoFCellAccessor::Implementation::
    set_dof_indices (*this, local_dof_indices);
}




template <class DH>
typename internal::DoFHandler::Iterators<DH>::cell_iterator
DoFCellAccessor<DH>::neighbor_child_on_subface (const unsigned int face,
						const unsigned int subface) const
{
  const TriaIterator<CellAccessor<dim,spacedim> > q
    = CellAccessor<dim,spacedim>::neighbor_child_on_subface (face, subface);
  return
    typename internal::DoFHandler::Iterators<DH>::cell_iterator (this->tria,
								 q->level (),
								 q->index (),
								 this->dof_handler);
}



template <class DH>
template <class InputVector, typename number>
void
DoFCellAccessor<DH>::
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
				     // if this cell has no children: simply
				     // return the exact values on this cell
    this->get_dof_values (values, interpolated_values);
  else
				     // otherwise clobber them from the children
    {
      Vector<number> tmp1(dofs_per_cell);
      Vector<number> tmp2(dofs_per_cell);

      interpolated_values = 0;

                                       // later on we will have to
                                       // push the values interpolated
                                       // from the child to the mother
                                       // cell into the output
                                       // vector. unfortunately, there
                                       // are two types of elements:
                                       // ones where you add up the
                                       // contributions from the
                                       // different child cells, and
                                       // ones where you overwrite.
                                       //
                                       // an example for the first is
                                       // piecewise constant (and
                                       // discontinuous) elements,
                                       // where we build the value on
                                       // the coarse cell by averaging
                                       // the values from the cell
                                       // (i.e. by adding up a
                                       // fraction of the values of
                                       // their values)
                                       //
                                       // an example for the latter
                                       // are the usual continuous
                                       // elements. the value on a
                                       // vertex of a coarse cell must
                                       // there be the same,
                                       // irrespective of the adjacent
                                       // cell we are presently on. so
                                       // we always overwrite. in
                                       // fact, we must, since we
                                       // cannot know in advance how
                                       // many neighbors there will
                                       // be, so there is no way to
                                       // compute the average with
                                       // fixed factors
                                       //
                                       // so we have to find out to
                                       // which type this element
                                       // belongs. the difficulty is:
                                       // the finite element may be a
                                       // composed one, so we can only
                                       // hope to do this for each
                                       // shape function
                                       // individually. in fact, there
                                       // are even weird finite
                                       // elements (for example the
                                       // Raviart-Thomas element)
                                       // which have shape functions
                                       // that are additive (interior
                                       // ones) and others that are
                                       // overwriting (face degrees of
                                       // freedom that need to be
                                       // continuous across the
                                       // face). to avoid checking
                                       // this over and over again, we
                                       // do this once now and cache
                                       // the results
      std::vector<bool> restriction_is_additive (dofs_per_cell);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        restriction_is_additive[i] = fe.restriction_is_additive(i);

      for (unsigned int child=0; child<this->n_children(); ++child)
	{
					   // get the values from the present
					   // child, if necessary by
					   // interpolation itself
	  this->child(child)->get_interpolated_dof_values (values,
							   tmp1);
					   // interpolate these to the mother
					   // cell
	  fe.get_restriction_matrix(child, this->refinement_case()).vmult (tmp2, tmp1);

                                           // and add up or set them
                                           // in the output vector
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
            if (restriction_is_additive[i])
              interpolated_values(i) += tmp2(i);
            else
              if (tmp2(i) != number())
                interpolated_values(i) = tmp2(i);
	}
    }
}



template <class DH>
template <class OutputVector, typename number>
void
DoFCellAccessor<DH>::
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
                                     // if this cell has no children: simply
				     // set the values on this cell
    this->set_dof_values (local_values, values);
  else
				     // otherwise distribute them to the children
    {
      Vector<number> tmp(dofs_per_cell);

      for (unsigned int child=0; child<this->n_children(); ++child)
	{
					   // prolong the given data
					   // to the present cell
	  this->get_fe().get_prolongation_matrix(child, this->refinement_case())
            .vmult (tmp, local_values);

	  this->child(child)->set_dof_values_by_interpolation (tmp, values);
	}
    }
}



// --------------------------------------------------------------------------
// explicit instantiations (for DoFHandler)

#if deal_II_dimension == 1
template class DoFAccessor<1, DoFHandler<1> >;
#endif

#if deal_II_dimension == 2
template class DoFAccessor<1, DoFHandler<2> >;
template class DoFAccessor<2, DoFHandler<2> >;

template class TriaRawIterator   <DoFAccessor<1, DoFHandler<2> > >;
template class TriaIterator      <DoFAccessor<1, DoFHandler<2> > >;
template class TriaActiveIterator<DoFAccessor<1, DoFHandler<2> > >;
#endif

#if deal_II_dimension == 3
template class DoFAccessor<1, DoFHandler<3> >;
template class DoFAccessor<2, DoFHandler<3> >;
template class DoFAccessor<3, DoFHandler<3> >;

template class TriaRawIterator   <DoFAccessor<1, DoFHandler<3> > >;
template class TriaIterator      <DoFAccessor<1, DoFHandler<3> > >;
template class TriaActiveIterator<DoFAccessor<1, DoFHandler<3> > >;
template class TriaRawIterator   <DoFAccessor<2, DoFHandler<3> > >;
template class TriaIterator      <DoFAccessor<2, DoFHandler<3> > >;
template class TriaActiveIterator<DoFAccessor<2, DoFHandler<3> > >;
#endif


template class DoFCellAccessor<DoFHandler<deal_II_dimension> >;

template class TriaRawIterator   <DoFCellAccessor<DoFHandler<deal_II_dimension> > >;
template class TriaIterator      <DoFCellAccessor<DoFHandler<deal_II_dimension> > >;
template class TriaActiveIterator<DoFCellAccessor<DoFHandler<deal_II_dimension> > >;


// --------------------------------------------------------------------------
// explicit instantiations (for hp::DoFHandler)

#if deal_II_dimension == 1
template class DoFAccessor<1, hp::DoFHandler<1> >;
#endif

#if deal_II_dimension == 2
template class DoFAccessor<1, hp::DoFHandler<2> >;
template class DoFAccessor<2, hp::DoFHandler<2> >;

template class TriaRawIterator   <DoFAccessor<1, hp::DoFHandler<2> > >;
template class TriaIterator      <DoFAccessor<1, hp::DoFHandler<2> > >;
template class TriaActiveIterator<DoFAccessor<1, hp::DoFHandler<2> > >;
#endif


#if deal_II_dimension == 3
template class DoFAccessor<1, hp::DoFHandler<3> >;
template class DoFAccessor<2, hp::DoFHandler<3> >;
template class DoFAccessor<3, hp::DoFHandler<3> >;

template class TriaRawIterator   <DoFAccessor<1, hp::DoFHandler<3> > >;
template class TriaIterator      <DoFAccessor<1, hp::DoFHandler<3> > >;
template class TriaActiveIterator<DoFAccessor<1, hp::DoFHandler<3> > >;
template class TriaRawIterator   <DoFAccessor<2, hp::DoFHandler<3> > >;
template class TriaIterator      <DoFAccessor<2, hp::DoFHandler<3> > >;
template class TriaActiveIterator<DoFAccessor<2, hp::DoFHandler<3> > >;
#endif


template class DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >;

template class TriaRawIterator   <DoFCellAccessor<hp::DoFHandler<deal_II_dimension> > >;
template class TriaIterator      <DoFCellAccessor<hp::DoFHandler<deal_II_dimension> > >;
template class TriaActiveIterator<DoFCellAccessor<hp::DoFHandler<deal_II_dimension> > >;



// // --------------------------------------------------------------------------
// // explicit instantiations (for DoFHandler)

#if deal_II_dimension == 1
template class DoFAccessor<1, DoFHandler<1,2> >;
#endif

#if deal_II_dimension == 2
template class DoFAccessor<1, DoFHandler<2,3> >;
template class DoFAccessor<2, DoFHandler<2,3> >;

template class TriaRawIterator   <DoFAccessor<1, DoFHandler<2,3> > >;
template class TriaIterator      <DoFAccessor<1, DoFHandler<2,3> > >;
template class TriaActiveIterator<DoFAccessor<1, DoFHandler<2,3> > >;
#endif


#if deal_II_dimension != 3
template class DoFCellAccessor<DoFHandler<deal_II_dimension,deal_II_dimension+1> >;

template class
TriaRawIterator   <DoFCellAccessor<DoFHandler<deal_II_dimension,deal_II_dimension+1> > >;
template class
TriaIterator      <DoFCellAccessor<DoFHandler<deal_II_dimension,deal_II_dimension+1> > >;
template class
TriaActiveIterator<DoFCellAccessor<DoFHandler<deal_II_dimension,deal_II_dimension+1> > >;
#endif

// --------------------------------------------------------------------------
// explicit instantiations (for hp::DoFHandler)

#if deal_II_dimension == 1
template class DoFAccessor<1, hp::DoFHandler<1,2> >;
#endif

#if deal_II_dimension == 2
template class DoFAccessor<1, hp::DoFHandler<2,3> >;
template class DoFAccessor<2, hp::DoFHandler<2,3> >;

template class TriaRawIterator   <DoFAccessor<1, hp::DoFHandler<2,3> > >;
template class TriaIterator      <DoFAccessor<1, hp::DoFHandler<2,3> > >;
template class TriaActiveIterator<DoFAccessor<1, hp::DoFHandler<2,3> > >;
#endif

#if deal_II_dimension != 3
template class DoFCellAccessor<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1> >;

template class
TriaRawIterator   <DoFCellAccessor<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1> > >;
template class
TriaIterator      <DoFCellAccessor<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1> > >;
template class
TriaActiveIterator<DoFCellAccessor<hp::DoFHandler<deal_II_dimension,deal_II_dimension+1> > >;
#endif


#include "dof_accessor.inst"

DEAL_II_NAMESPACE_CLOSE
