//----------------------------  dof_iterator_selector.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_iterator_selector.h  ---------------------------
#ifndef __deal2__dof_iterator_selector_h
#define __deal2__dof_iterator_selector_h

template <class DH> class DoFAccessor;
template <class DH> class DoFCellAccessor;
template <int celldim, class DH> class DoFObjectAccessor;
template <int dim> class FiniteElement;
template <int dim, typename Accessor> class TriaRawIterator;
template <int dim, typename Accessor> class TriaIterator;
template <int dim, typename Accessor> class TriaActiveIterator;
template <int dim> class Triangulation;
template <int dim> class DoFHandler;


namespace internal 
{
/**
 * Define some types which differ between the dimensions. This class
 * is analogous to the TriaDimensionInfo class hierarchy.
 * 
 * @ref DoFIteratorSelector<1>
 * @ref DoFIteratorSelector<2>
 * @ref DoFIteratorSelector<3>
 *
 * @author Wolfgang Bangerth, 1998; Oliver Kayser-Herold and Wolfgang Bangerth, 2003
 */
  template <class DH>
  struct DoFIteratorSelector
  {
  };



/**
 * Define some types for the DoF handling in one dimension.
 *
 * The types have the same meaning as those declared in
 * TriaDimensionInfo<2>, only the treatment of templates is a
 * little more complicated since the @p DoFAccessor classes want a
 * template template argument.
 *
 * @author Wolfgang Bangerth, 1998; Oliver Kayser-Herold and Wolfgang Bangerth, 2003
 */
  template <template <int> class DH>
  struct DoFIteratorSelector<DH<1> >
  {
      typedef TriaRawIterator<1,DoFCellAccessor<DH<1> > >    raw_line_iterator;
      typedef TriaIterator<1,DoFCellAccessor<DH<1> > >       line_iterator;
      typedef TriaActiveIterator<1,DoFCellAccessor<DH<1> > > active_line_iterator;
      
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
 * Define some types for the DoF handling in two dimensions.
 *
 * The types have the same meaning as those declared in
 * TriaDimensionInfo<2>, only the treatment of templates is a
 * little more complicated since the @p DoFAccessor classes want a
 * template template argument.
 *
 * @author Wolfgang Bangerth, 1998; Oliver Kayser-Herold and Wolfgang Bangerth, 2003
 */
  template <template <int> class DH>
  struct DoFIteratorSelector<DH<2> >
  {
      typedef TriaRawIterator<2,DoFObjectAccessor<1, DH<2> > >    raw_line_iterator;
      typedef TriaIterator<2,DoFObjectAccessor<1, DH<2> > >       line_iterator;
      typedef TriaActiveIterator<2,DoFObjectAccessor<1, DH<2> > > active_line_iterator;
      
      typedef TriaRawIterator<2,DoFCellAccessor<DH<2> > >         raw_quad_iterator;
      typedef TriaIterator<2,DoFCellAccessor<DH<2> > >            quad_iterator;
      typedef TriaActiveIterator<2,DoFCellAccessor<DH<2> > >      active_quad_iterator;
      
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
 * Define some types for the DoF handling in two dimensions.
 *
 * The types have the same meaning as those declared in
 * TriaDimensionInfo<3>, only the treatment of templates is a
 * little more complicated since the @p DoFAccessor classes want a
 * template template argument.
 *
 * @author Wolfgang Bangerth, 1998; Oliver Kayser-Herold and Wolfgang Bangerth, 2003
 */
  template <template <int> class DH>
  struct DoFIteratorSelector<DH<3> >
  {
      typedef TriaRawIterator<3,DoFObjectAccessor<1, DH<3> > >    raw_line_iterator;
      typedef TriaIterator<3,DoFObjectAccessor<1, DH<3> > >       line_iterator;
      typedef TriaActiveIterator<3,DoFObjectAccessor<1, DH<3> > > active_line_iterator;
      
      typedef TriaRawIterator<3,DoFObjectAccessor<2, DH<3> > >    raw_quad_iterator;
      typedef TriaIterator<3,DoFObjectAccessor<2, DH<3> > >       quad_iterator;
      typedef TriaActiveIterator<3,DoFObjectAccessor<2, DH<3> > > active_quad_iterator;
      
      typedef TriaRawIterator<3,DoFCellAccessor<DH<3> > >         raw_hex_iterator;
      typedef TriaIterator<3,DoFCellAccessor<DH<3> > >            hex_iterator;
      typedef TriaActiveIterator<3,DoFCellAccessor<DH<3> > >      active_hex_iterator;
      
      typedef raw_hex_iterator    raw_cell_iterator;
      typedef hex_iterator        cell_iterator;
      typedef active_hex_iterator active_cell_iterator;
      
      typedef raw_quad_iterator    raw_face_iterator;
      typedef quad_iterator        face_iterator;
      typedef active_quad_iterator active_face_iterator;    
  };
}

#endif // __deal2__dof_iterator_selector_h
