//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_dof_iterator_selector_h
#define __deal2__mg_dof_iterator_selector_h


#include <base/config.h>

template <int celldim, int dim> class MGDoFAccessor;
template <int celldim, int dim> class MGDoFObjectAccessor;
template <int dim>              class MGDoFObjectAccessor<0, dim>;
template <int dim>              class MGDoFObjectAccessor<1, dim>;
template <int dim>              class MGDoFObjectAccessor<2, dim>;
template <int dim>              class MGDoFObjectAccessor<3, dim>;
template <int dim> class MGDoFCellAccessor;

template <int dim, typename Accessor> class TriaRawIterator;
template <int dim, typename Accessor> class TriaIterator;
template <int dim, typename Accessor> class TriaActiveIterator;


namespace internal
{
  namespace MGDoFHandler
  {
/*!@addtogroup mg */
/*@{*/

/**
 * Depending on the dimension, this class defines the iterator types for
 * MGDoFHandler objects.
 *
 * @note The typedefs in this class are only for reference. Their real
 * values are to be found in the specializations for different
 * dimensions.
 */
    template <int dim>
    class Iterators
    {
                                         /// Line iterator for internal use
        typedef void * raw_line_iterator;
                                         /// Line iterator for all lines of hierarchy
        typedef void * line_iterator;
                                         /// Iterator for active lines
        typedef void * active_line_iterator;
    
                                         /// Quad iterator for internal use
        typedef void * raw_quad_iterator;
                                         /// Quad iterator for all quadrilaterals of hierarchy
        typedef void * quad_iterator;
                                         /// Iterator for active quadrilaterals
        typedef void * active_quad_iterator;

                                         /// Hex iterator for internal use
        typedef void * raw_hex_iterator;
                                         /// Hex iterator for all hexahedra of hierarchy
        typedef void * hex_iterator;
                                         /// Iterator for active hexahedra
        typedef void * active_hex_iterator;

                                         /** Cell iterator for internal
                                          * use. It is defined as one of
                                          * the raw iterators above,
                                          * depending on the dimension.
                                          */
        typedef void * raw_cell_iterator;
                                         /** Cell iterator for all cells
                                          * of hierarchy. It is defined as
                                          * one of the iterators above,
                                          * depending on the dimension.
                                          */
        typedef void * cell_iterator;
                                         /** Iterator for active cells. It
                                          * is defined as one of the
                                          * active iterators above,
                                          * depending on the dimension.
                                          */
        typedef void * active_cell_iterator;

                                         /** Face iterator for internal
                                          * use. It is defined as one of
                                          * the raw iterators above,
                                          * depending on the dimension.
                                          */
        typedef void * raw_face_iterator;
                                         /** Face iterator for all faces
                                          * of hierarchy. It is defined as
                                          * one of the iterators above,
                                          * depending on the dimension.
                                          */
        typedef void * face_iterator;
                                         /** Iterator for active faces. It
                                          * is defined as one of the
                                          * active iterators above,
                                          * depending on the dimension.
                                          */
        typedef void * active_face_iterator;
    };



/**
 * Iterators for MGDofHandler in one dimension. See the template
 * Iterators for further explanation.
 */
    template <>
    class Iterators<1>
    {
      public:
        typedef TriaRawIterator<1,MGDoFCellAccessor<1> >    raw_line_iterator;
        typedef TriaIterator<1,MGDoFCellAccessor<1> >       line_iterator;
        typedef TriaActiveIterator<1,MGDoFCellAccessor<1> > active_line_iterator;

        typedef void * raw_quad_iterator;
        typedef void * quad_iterator;
        typedef void * active_quad_iterator;

        typedef void * raw_hex_iterator;
        typedef void * hex_iterator;
        typedef void * active_hex_iterator;

        typedef raw_line_iterator    raw_cell_iterator;
        typedef line_iterator        cell_iterator;
        typedef active_line_iterator active_cell_iterator;

        typedef void * raw_face_iterator;
        typedef void * face_iterator;
        typedef void * active_face_iterator;
    };



/**
 * Iterators for MGDofHandler in one dimension. See the template
 * Iterators for further explanation.
 */
    template <>
    class Iterators<2>
    {
      public:
        typedef TriaRawIterator<2,MGDoFObjectAccessor<1, 2> >    raw_line_iterator;
        typedef TriaIterator<2,MGDoFObjectAccessor<1, 2> >       line_iterator;
        typedef TriaActiveIterator<2,MGDoFObjectAccessor<1, 2> > active_line_iterator;
    
        typedef TriaRawIterator<2,MGDoFCellAccessor<2> >         raw_quad_iterator;
        typedef TriaIterator<2,MGDoFCellAccessor<2> >            quad_iterator;
        typedef TriaActiveIterator<2,MGDoFCellAccessor<2> >      active_quad_iterator;

        typedef void * raw_hex_iterator;
        typedef void * hex_iterator;
        typedef void * active_hex_iterator;

        typedef raw_quad_iterator    raw_cell_iterator;
        typedef quad_iterator        cell_iterator;
        typedef active_quad_iterator active_cell_iterator;

        typedef raw_line_iterator    raw_face_iterator;
        typedef line_iterator        face_iterator;
        typedef active_line_iterator active_face_iterator;    
    };



/**
 * Iterators for MGDofHandler in one dimension. See the template
 * Iterators for further explanation.
 */
    template <>
    class Iterators<3>
    {
      public:
        typedef TriaRawIterator<3,MGDoFObjectAccessor<1, 3> >    raw_line_iterator;
        typedef TriaIterator<3,MGDoFObjectAccessor<1, 3> >       line_iterator;
        typedef TriaActiveIterator<3,MGDoFObjectAccessor<1, 3> > active_line_iterator;

        typedef TriaRawIterator<3,MGDoFObjectAccessor<2, 3> >    raw_quad_iterator;
        typedef TriaIterator<3,MGDoFObjectAccessor<2, 3> >       quad_iterator;
        typedef TriaActiveIterator<3,MGDoFObjectAccessor<2, 3> > active_quad_iterator;

        typedef TriaRawIterator<3,MGDoFCellAccessor<3> >               raw_hex_iterator;
        typedef TriaIterator<3,MGDoFCellAccessor<3> >                  hex_iterator;
        typedef TriaActiveIterator<3,MGDoFCellAccessor<3> >            active_hex_iterator;

        typedef raw_hex_iterator    raw_cell_iterator;
        typedef hex_iterator        cell_iterator;
        typedef active_hex_iterator active_cell_iterator;

        typedef raw_quad_iterator    raw_face_iterator;
        typedef quad_iterator        face_iterator;
        typedef active_quad_iterator active_face_iterator;    
    };
  }
}


#endif
/*----------------------------   mg_dof_iterator_selector.h     ---------------------------*/
