//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_dof_iterator_selector_h
#define __deal2__mg_dof_iterator_selector_h


#include <deal.II/base/config.h>


DEAL_II_NAMESPACE_OPEN

template <int structdim, int dim, int spacedim> class MGDoFAccessor;
template <int dim, int spacedim> class MGDoFCellAccessor;
template <int, int, int> class InvalidAccessor;

template <typename Accessor> class TriaRawIterator;
template <typename Accessor> class TriaIterator;
template <typename Accessor> class TriaActiveIterator;


namespace internal
{
  namespace MGDoFHandler
  {
    template <int dim, int spacedim>
    class Iterators;


/**
 * Iterators for MGDofHandler in one dimension. See the @ref Iterators module
 * for more information.
 */
    template <int spacedim>
    class Iterators<1,spacedim>
    {
      public:
	typedef dealii::MGDoFCellAccessor<1,spacedim> CellAccessor;
	typedef dealii::MGDoFAccessor<0,1,spacedim> FaceAccessor;

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
 * Iterators for MGDofHandler in two dimensions. See the @ref Iterators module
 * for more information.
 */
    template <int spacedim>
    class Iterators<2,spacedim>
    {
      public:
	typedef dealii::MGDoFCellAccessor<2,spacedim> CellAccessor;
	typedef dealii::MGDoFAccessor<1,2,spacedim> FaceAccessor;

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
 * Iterators for MGDofHandler in three dimensions. See the @ref Iterators module
 * for more information.
 */
    template <int spacedim>
    class Iterators<3,spacedim>
    {
      public:
	typedef dealii::MGDoFCellAccessor<3,spacedim> CellAccessor;
	typedef dealii::MGDoFAccessor<2,3,spacedim> FaceAccessor;

        typedef TriaRawIterator   <dealii::MGDoFAccessor<1,3,spacedim> >    raw_line_iterator;
        typedef TriaIterator      <dealii::MGDoFAccessor<1,3,spacedim> >       line_iterator;
        typedef TriaActiveIterator<dealii::MGDoFAccessor<1,3,spacedim> > active_line_iterator;

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

#endif
/*----------------------------   mg_dof_iterator_selector.h     ---------------------------*/
