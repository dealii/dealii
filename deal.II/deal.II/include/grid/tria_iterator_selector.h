//----------------------------  tria_iterator_selector.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tria_iterator_selector.h  ---------------------------
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

/**
 *  This class implements some types which differ between the
 *  dimensions.  Declare it to have a template parameter, but do not
 *  actually declare anything concrete apart from the other classes
 *  which are explicitly instantiated ones with the same name.
 *
 * @author Wolfgang Bangerth, 1998
 */
  template <int dim>
  struct TriaIteratorSelector
  {
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
  struct TriaIteratorSelector<1>
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
  struct TriaIteratorSelector<2>
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
 *  as for lower dimensions (see <tt>TriaIteratorSelector<[12]></tt>). The
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
  struct TriaIteratorSelector<3>
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

#endif // __deal2__tria_iterator_selector_h
