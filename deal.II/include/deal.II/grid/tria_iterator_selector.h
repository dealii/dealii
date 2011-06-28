//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tria_iterator_selector_h
#define __deal2__tria_iterator_selector_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class CellAccessor;
template <int, int, int> class TriaAccessorBase;
template <int, int, int> class InvalidAccessor;
template <int, int, int> class TriaAccessor;
template <int dim, int spacedim>  class TriaAccessor<0, dim, spacedim>;
template <typename Accessor> class TriaRawIterator;
template <typename Accessor> class TriaIterator;
template <typename Accessor> class TriaActiveIterator;

namespace internal 
{
  namespace Triangulation
  {
    template <int dim, int spacedim>
    struct Iterators;

/**
 *  This class implements some types which differ between the dimensions.
 *  These are the declararions for the 1D case only. See the @ref Iterators
 *  module for more information.
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
 *  iterators over InvalidAccessor. Thus these types exist, but
 *  are useless and will certainly make any involuntary use
 *  visible. The same holds for hexahedron iterators.
 *
 *  The same applies for the @p face_iterator types, since lines
 *  have no substructures apart from vertices, which are handled in
 *  a different way, however.
 *
 * @author Wolfgang Bangerth, 1998
 */
    template <int spacedim>
    struct Iterators<1,spacedim>
    {
	typedef dealii::CellAccessor<1,spacedim> CellAccessor;
	typedef InvalidAccessor<0,1,spacedim> FaceAccessor;
	
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
 *  This class implements some types which differ between the dimensions.
 *  These are the declararions for the 2D case only. See the @ref Iterators
 *  module for more information.
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
 *  iterators over InvalidAccessor. Thus these types exist, but
 *  are useless and will certainly make any involuntary use visible.
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
    template <int spacedim>
    struct Iterators<2,spacedim>
    {
	typedef dealii::CellAccessor<2, spacedim> CellAccessor;
	typedef dealii::TriaAccessor<1, 2, spacedim> FaceAccessor;
	
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
 *  This class implements some types which differ between the dimensions.
 *  These are the declararions for the 3D case only. See the @ref Iterators
 *  module for more information.
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
    template <int spacedim>
    struct Iterators<3,spacedim>
    {
	typedef dealii::CellAccessor<3, spacedim> CellAccessor;
	typedef dealii::TriaAccessor<2, 3, spacedim> FaceAccessor;
	
        typedef TriaRawIterator   <dealii::TriaAccessor<1, 3, spacedim> > raw_line_iterator;
        typedef TriaIterator      <dealii::TriaAccessor<1, 3, spacedim> > line_iterator;
        typedef TriaActiveIterator<dealii::TriaAccessor<1, 3, spacedim> > active_line_iterator;
    
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

#endif // __deal2__tria_iterator_selector_h
