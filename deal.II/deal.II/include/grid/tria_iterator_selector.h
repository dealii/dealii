//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tria_iterator_selector_h
#define __deal2__tria_iterator_selector_h


template <int dim> class CellAccessor;
template <int dim> class TriaAccessor;
template <int, int> class TriaObjectAccessor;
template <int dim>  class TriaObjectAccessor<0, dim>;
template <int dim>  class TriaObjectAccessor<1, dim>;
template <int dim>  class TriaObjectAccessor<2, dim>;
template <int dim>  class TriaObjectAccessor<3, dim>;
template <int dim, typename Accessor> class TriaRawIterator;
template <int dim, typename Accessor> class TriaIterator;
template <int dim, typename Accessor> class TriaActiveIterator;


namespace internal 
{
  namespace Triangulation
  {
/**
 * Definition of the iterator types of the Triangulation.
 *
 * @note The actual definitions used are defined in specializations of
 * this class. The template is just here for documentation and shows a
 * generic case, while it is clear that for instance #quad_iterator is
 * of no use in one dimension.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 2006
 */
    template <int dim>
    struct Iterators
    {
                                         /// Iterate on raw lines
        typedef TriaRawIterator<dim,TriaObjectAccessor<1, dim> > raw_line_iterator;
                                         /// Iterate on raw lines
        typedef TriaIterator<dim,TriaObjectAccessor<1, dim> > line_iterator;
                                         /// Iterate on raw lines
        typedef TriaActiveIterator<dim,TriaObjectAccessor<1, dim> > active_line_iterator;
    
                                         /// Iterate on raw quadrilaterals
        typedef TriaRawIterator<dim,TriaObjectAccessor<2, dim> > raw_quad_iterator;
                                         /// Iterate on quadrilaterals
        typedef TriaIterator<dim,TriaObjectAccessor<2, dim> > quad_iterator;
                                         /// Iterate on active quadrilaterals
        typedef TriaActiveIterator<dim,TriaObjectAccessor<2, dim> > active_quad_iterator;

                                         /// Iterate on raw hexahedra
        typedef TriaRawIterator<dim,TriaObjectAccessor<3, dim> > raw_hex_iterator;
                                         /// Iterate on hexahedra
        typedef TriaIterator<dim,TriaObjectAccessor<3, dim> > hex_iterator;
                                         /// Iterate on active hexahedra
        typedef TriaActiveIterator<dim,TriaObjectAccessor<3, dim> > active_hex_iterator;

                                         /// Iterate on raw cells
        typedef TriaRawIterator<dim,CellAccessor<dim> > raw_cell_iterator;
                                         /// Iterate on cells
        typedef TriaIterator<dim,CellAccessor<dim> > cell_iterator;
                                         /// Iterate on active cells
        typedef TriaActiveIterator<dim,CellAccessor<dim> > active_cell_iterator;

                                         /// Iterate on raw faces
        typedef TriaRawIterator<dim,TriaObjectAccessor<dim-1, dim> > raw_face_iterator;
                                         /// Iterate on faces
        typedef TriaIterator<dim,TriaObjectAccessor<dim-1, dim> > face_iterator;
                                         /// Iterate on active faces
        typedef TriaActiveIterator<dim,TriaObjectAccessor<dim-1, dim> > active_face_iterator;
    };



/**
 *  This class implements some types which differ between the dimensions.
 *  These are the declararions for the 1D case only.
 *
 *  A @p line_iterator is typdef'd to an iterator operating on the
 *  @p lines member variable of a <tt>Triangulation<1></tt> object. An
 *  @p active_line_iterator only operates on the active lines.
 *  @p raw_line_iterator objects operate on all lines, used or not.
 *
 *  Since we are in one dimension, the following identities are declared:
 *  @verbatim
 *    typedef raw_line_iterator    raw_cell_iterator;
 *    typedef line_iterator        cell_iterator;
 *    typedef active_line_iterator active_cell_iterator;
 *  @endverbatim
 *
 *  To enable the declaration of @p begin_quad and the like in
 *  <tt>Triangulation<1></tt>, the @p quad_iterators are declared as
 *  <tt>void *</tt>. Thus these types exist, but are useless and will
 *  certainly make any involuntary use visible. The same holds
 *  for hexahedron iterators.
 *
 *  The same applies for the @p face_iterator types, since lines
 *  have no substructures apart from vertices, which are handled in
 *  a different way, however.
 *
 * @author Wolfgang Bangerth, 1998
 */
    template <>
    struct Iterators<1>
    {
        typedef TriaRawIterator<1,CellAccessor<1> >    raw_line_iterator;
        typedef TriaIterator<1,CellAccessor<1> >       line_iterator;
        typedef TriaActiveIterator<1,CellAccessor<1> > active_line_iterator;

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
 *  This class implements some types which differ between the dimensions.
 *  These are the declararions for the 2D case only.
 *
 *  A @p line_iterator is typdef'd to an iterator operating on the
 *  @p lines member variable of a <tt>Triangulation<2></tt> object. An
 *  @p active_line_iterator only operates on the active lines.
 *  @p raw_line_iterator objects operate on all lines, used or not.
 *  Using @p active_line_iterators may not be particularly in 2D useful since it
 *  only operates on unrefined lines. However, also refined lines may bound
 *  unrefined cells if the neighboring cell is refined once more than the
 *  present one.
 *
 *  Similarly to line iterators, @p quad_iterator, @p raw_quad_iterator and
 *  @p active_quad_iterator are declared.
 *  
 *  To enable the declaration of @p begin_hex and the like in
 *  <tt>Triangulation<[12]></tt>, the @p hex_iterators are declared as
 *  <tt>void *</tt>. Thus these types exist, but are useless and will
 *  certainly make any involuntary use visible.
 *
 *  Since we are in two dimension, the following identities are declared:
 *  @verbatim
 *    typedef raw_quad_iterator    raw_cell_iterator;
 *    typedef quad_iterator        cell_iterator;
 *    typedef active_quad_iterator active_cell_iterator;
 *
 *    typedef raw_line_iterator    raw_face_iterator;
 *    typedef line_iterator        face_iterator;
 *    typedef active_line_iterator active_face_iterator;    
 *  @endverbatim
 *
 * @author Wolfgang Bangerth, 1998
 */
    template <>
    struct Iterators<2>
    {
        typedef TriaRawIterator<2,TriaObjectAccessor<1, 2> >    raw_line_iterator;
        typedef TriaIterator<2,TriaObjectAccessor<1, 2> >       line_iterator;
        typedef TriaActiveIterator<2,TriaObjectAccessor<1, 2> > active_line_iterator;
    
        typedef TriaRawIterator<2,CellAccessor<2> >    raw_quad_iterator;
        typedef TriaIterator<2,CellAccessor<2> >       quad_iterator;
        typedef TriaActiveIterator<2,CellAccessor<2> > active_quad_iterator;

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
 *  This class implements some types which differ between the dimensions.
 *  These are the declararions for the 3D case only.
 *
 *  For the declarations of the data types, more or less the same holds
 *  as for lower dimensions (see <tt>Iterators<[12]></tt>). The
 *  dimension specific data types are here, since we are in three dimensions:
 *  @verbatim
 *    typedef raw_hex_iterator    raw_cell_iterator;
 *    typedef hex_iterator        cell_iterator;
 *    typedef active_hex_iterator active_cell_iterator;
 *
 *    typedef raw_quad_iterator    raw_face_iterator;
 *    typedef quad_iterator        face_iterator;
 *    typedef active_quad_iterator active_face_iterator;    
 *  @endverbatim
 *
 * @author Wolfgang Bangerth, 1998
 */
    template <>
    struct Iterators<3>
    {
        typedef TriaRawIterator<3,TriaObjectAccessor<1, 3> >    raw_line_iterator;
        typedef TriaIterator<3,TriaObjectAccessor<1, 3> >       line_iterator;
        typedef TriaActiveIterator<3,TriaObjectAccessor<1, 3> > active_line_iterator;
    
        typedef TriaRawIterator<3,TriaObjectAccessor<2, 3> >    raw_quad_iterator;
        typedef TriaIterator<3,TriaObjectAccessor<2, 3> >       quad_iterator;
        typedef TriaActiveIterator<3,TriaObjectAccessor<2, 3> > active_quad_iterator;

        typedef TriaRawIterator<3,CellAccessor<3> >    raw_hex_iterator;
        typedef TriaIterator<3,CellAccessor<3> >       hex_iterator;
        typedef TriaActiveIterator<3,CellAccessor<3> > active_hex_iterator;

        typedef raw_hex_iterator    raw_cell_iterator;
        typedef hex_iterator        cell_iterator;
        typedef active_hex_iterator active_cell_iterator;

        typedef raw_quad_iterator    raw_face_iterator;
        typedef quad_iterator        face_iterator;
        typedef active_quad_iterator active_face_iterator;    
    };
  
  }
  
}

#endif // __deal2__tria_iterator_selector_h
