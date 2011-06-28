//----------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2006, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------
#ifndef __deal2__dof_iterators_h
#define __deal2__dof_iterators_h

#include <deal.II/base/config.h>


DEAL_II_NAMESPACE_OPEN

template <int, int, int> class InvalidAccessor;

template <int structdim, class DH> class DoFAccessor;
template <class DH> class DoFCellAccessor;

template <int dim, int spacedim> class FiniteElement;
template <typename Accessor> class TriaRawIterator;
template <typename Accessor> class TriaIterator;
template <typename Accessor> class TriaActiveIterator;
template <int dim, int spacedim> class Triangulation;
template <int dim, int spacedim> class DoFHandler;


namespace internal
{
  namespace DoFHandler
  {
    template <class DH>
    struct Iterators;


/**
 * Define some types for DoF handling in one dimension.
 *
 * The types have the same meaning as those declared in
 * internal::Triangulation::Iterators<1,spacedim>, only the treatment of
 * templates is a little more complicated. See the @ref Iterators module
 * for more information.
 *
 * @author Wolfgang Bangerth, Oliver Kayser-Herold, Guido Kanschat, 1998, 2003, 2008, 2010
 */
    template <template <int, int> class DH, int spacedim>
    struct Iterators<DH<1, spacedim> >
    {
        typedef DH<1,spacedim> DoFHandler_type;
	typedef dealii::DoFCellAccessor<DoFHandler_type> CellAccessor;
	typedef dealii::DoFAccessor<0,DoFHandler_type> FaceAccessor;

        typedef TriaRawIterator   <CellAccessor> raw_line_iterator;
        typedef TriaIterator      <CellAccessor> line_iterator;
        typedef TriaActiveIterator<CellAccessor> active_line_iterator;

        typedef TriaRawIterator   <InvalidAccessor<2,1,spacedim> > raw_quad_iterator;
        typedef TriaIterator      <InvalidAccessor<2,1,spacedim> > quad_iterator;
        typedef TriaActiveIterator<InvalidAccessor<2,1,spacedim> > active_quad_iterator;

        typedef TriaRawIterator   <InvalidAccessor<3,1,spacedim> > raw_hex_iterator;
        typedef TriaIterator      <InvalidAccessor<3,1,spacedim> > hex_iterator;
        typedef TriaActiveIterator<InvalidAccessor<3,1,spacedim> > active_hex_iterator;

        typedef raw_line_iterator    raw_cell_iterator;
        typedef line_iterator        cell_iterator;
        typedef active_line_iterator active_cell_iterator;

        typedef TriaRawIterator   <FaceAccessor> raw_face_iterator;
        typedef TriaIterator      <FaceAccessor> face_iterator;
        typedef TriaActiveIterator<FaceAccessor> active_face_iterator;
    };




/**
 * Define some types for DoF handling in two dimensions.
 *
 * The types have the same meaning as those declared in
 * internal::Triangulation::Iterators<2,spacedim>, only the treatment of
 * templates is a little more complicated. See the @ref Iterators module
 * for more information.
 *
 * @author Wolfgang Bangerth, Oliver Kayser-Herold, Guido Kanschat, 1998, 2003, 2008, 2010
 */
    template <template <int, int> class DH, int spacedim>
    struct Iterators<DH<2, spacedim> >
    {
        typedef DH<2,spacedim> DoFHandler_type;
	typedef dealii::DoFCellAccessor<DoFHandler_type> CellAccessor;
	typedef dealii::DoFAccessor<1, DoFHandler_type> FaceAccessor;

        typedef TriaRawIterator   <FaceAccessor> raw_line_iterator;
        typedef TriaIterator      <FaceAccessor> line_iterator;
        typedef TriaActiveIterator<FaceAccessor> active_line_iterator;

        typedef TriaRawIterator   <CellAccessor> raw_quad_iterator;
        typedef TriaIterator      <CellAccessor> quad_iterator;
        typedef TriaActiveIterator<CellAccessor> active_quad_iterator;

        typedef TriaRawIterator   <InvalidAccessor<3,2,spacedim> > raw_hex_iterator;
        typedef TriaIterator      <InvalidAccessor<3,2,spacedim> > hex_iterator;
        typedef TriaActiveIterator<InvalidAccessor<3,2,spacedim> > active_hex_iterator;

        typedef raw_quad_iterator    raw_cell_iterator;
        typedef quad_iterator        cell_iterator;
        typedef active_quad_iterator active_cell_iterator;

        typedef raw_line_iterator    raw_face_iterator;
        typedef line_iterator        face_iterator;
        typedef active_line_iterator active_face_iterator;
    };




/**
 * Define some types for DoF handling in three dimensions.
 *
 * The types have the same meaning as those declared in
 * internal::Triangulation::Iterators<3,spacedim>, only the treatment of
 * templates is a little more complicated. See the @ref Iterators module
 * for more information.
 *
 * @author Wolfgang Bangerth, Oliver Kayser-Herold, Guido Kanschat, 1998, 2003, 2008, 2010
 */
    template <template <int, int> class DH, int spacedim>
    struct Iterators<DH<3, spacedim> >
    {
        typedef DH<3, spacedim> DoFHandler_type;
	typedef dealii::DoFCellAccessor<DoFHandler_type> CellAccessor;
	typedef dealii::DoFAccessor<2, DoFHandler_type> FaceAccessor;

        typedef TriaRawIterator   <dealii::DoFAccessor<1, DoFHandler_type> > raw_line_iterator;
        typedef TriaIterator      <dealii::DoFAccessor<1, DoFHandler_type> > line_iterator;
        typedef TriaActiveIterator<dealii::DoFAccessor<1, DoFHandler_type> > active_line_iterator;

        typedef TriaRawIterator   <FaceAccessor> raw_quad_iterator;
        typedef TriaIterator      <FaceAccessor> quad_iterator;
        typedef TriaActiveIterator<FaceAccessor> active_quad_iterator;

        typedef TriaRawIterator   <CellAccessor> raw_hex_iterator;
        typedef TriaIterator      <CellAccessor> hex_iterator;
        typedef TriaActiveIterator<CellAccessor> active_hex_iterator;

        typedef raw_hex_iterator    raw_cell_iterator;
        typedef hex_iterator        cell_iterator;
        typedef active_hex_iterator active_cell_iterator;

        typedef raw_quad_iterator    raw_face_iterator;
        typedef quad_iterator        face_iterator;
        typedef active_quad_iterator active_face_iterator;
    };
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif // __deal2__dof_iterator_selector_h
