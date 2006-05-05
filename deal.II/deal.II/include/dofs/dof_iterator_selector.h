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
#ifndef __deal2__dof_iterators_h
#define __deal2__dof_iterators_h

//TODO: Rename this file to dof_iterators.h
//TODO: Can we avoid the specializations?

// To consider for the second TODO:
//  1. Define accessors for unreasonable dimensions
//  2. Can cell_iterator and quad_iterator in 2D be different (or hex
//     in 3D)?

template <int structdim, class DH> class DoFAccessor;
template <class DH> class DoFCellAccessor;
template <int celldim, class DH> class DoFObjectAccessor;
template <int dim> class FiniteElement;
template <int dim, typename Accessor> class TriaRawIterator;
template <int dim, typename Accessor> class TriaIterator;
template <int dim, typename Accessor> class TriaActiveIterator;
template <int dim> class Triangulation;
template <int dim> class DoFHandler;


/**
 * A pseudo class defining the iterator types used by DoFHandler and
 * hp::DoFHandler.  The typedefs in this class are synonymous with
 * those in the handler classes, such that this class can be used in
 * order to avoid including the complete header files of the handler.
 *
 * The class template is not used, but only its specializations with
 * respect to the dimension. These differ from the template in two
 * points: first, the iterator to objects of the same dimension as the
 * handler is aliased to the #cell_iterator. Second, #face_iterator is
 * aliased to objects of codimension one.
 *
 *
 * @warning This raw iterators are not normally used in application
 * programs.
 *
 * @ingroup dofs
 * @ingroup Accessors
 * @author W. Bangerth, G. Kanschat,  O. Kayser-Herold, 1998, 2003, 2006
 */
template <class DH>
struct DoFIterators
{
				     /// The dof handler class.
    typedef DH DoFHandler_type;
				     /// The topological dimension of the dof handler.
    static const unsigned int dim = DH::dimension;
    
				     /// Iterator for raw lines.
    typedef TriaRawIterator<dim, DoFObjectAccessor<1, DH> > raw_line_iterator;
				     /// Iterator for usual lines.
    typedef TriaIterator<dim, DoFObjectAccessor<1, DH> > line_iterator;
				     /// Iterator for active lines.
    typedef TriaActiveIterator<dim, DoFObjectAccessor<1, DH> > active_line_iterator;
				     /// Iterator for raw quadrilaterals
    typedef TriaRawIterator<dim, DoFObjectAccessor<2, DH> > raw_quad_iterator;
				     /// Iterator for quadrilaterals
    typedef TriaIterator<dim ,DoFObjectAccessor<2, DH> > quad_iterator;
				     /// Iterator for active quadrilaterals
    typedef TriaActiveIterator<dim ,DoFObjectAccessor<2, DH> > active_quad_iterator;
				     /// Iterator for raw hexahedra
    typedef TriaRawIterator<dim ,DoFObjectAccessor<3, DH> > raw_hex_iterator;
				     /// Iterator for hexahedra
    typedef TriaIterator<dim ,DoFObjectAccessor<3, DH> > hex_iterator;
				     /// Iterator for active hexahedra
    typedef TriaActiveIterator<dim ,DoFObjectAccessor<3, DH> > active_hex_iterator;
				     /// Iterator for raw cells
    typedef TriaRawIterator<dim, DoFCellAccessor<DH> > raw_cell_iterator;
				     /// Iterator for cells
    typedef TriaIterator<dim, DoFCellAccessor<DH> > cell_iterator;
				     /// Iterator for active cells
    typedef TriaActiveIterator<dim, DoFCellAccessor<DH> > active_cell_iterator;
				     /// Iterator for raw faces
    typedef TriaRawIterator<dim ,DoFObjectAccessor<dim-1, DH> > raw_face_iterator;
				     /// Iterator for faces
    typedef TriaIterator<dim ,DoFObjectAccessor<dim-1, DH> > face_iterator;
				     /// Iterator for active faces
    typedef TriaActiveIterator<dim ,DoFObjectAccessor<dim-1, DH> > active_face_iterator;
};



/**
 * Define some types for the DoF handling in one dimension.
 *
 * The types have the same meaning as those declared in
 * TriaDimensionInfo<2>, only the treatment of templates is a
 * little more complicated since the @p DoFAccessor classes want a
 * template template argument.
 *
 * @author Wolfgang Bangerth, Oliver Kayser-Herold, 1998, 2003
 */
template <template <int> class DH>
struct DoFIterators<DH<1> >
{
    typedef DH<1> DoFHandler_type;
    
    typedef TriaRawIterator<1,DoFCellAccessor<DoFHandler_type> >    raw_line_iterator;
    typedef TriaIterator<1,DoFCellAccessor<DoFHandler_type> >       line_iterator;
    typedef TriaActiveIterator<1,DoFCellAccessor<DoFHandler_type> > active_line_iterator;
    
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
 * @author Wolfgang Bangerth, Oliver Kayser-Herold, 1998, 2003
 */
template <template <int> class DH>
struct DoFIterators<DH<2> >
{
    typedef DH<2> DoFHandler_type;
    
    typedef TriaRawIterator<2,DoFObjectAccessor<1, DoFHandler_type> >    raw_line_iterator;
    typedef TriaIterator<2,DoFObjectAccessor<1, DoFHandler_type> >       line_iterator;
    typedef TriaActiveIterator<2,DoFObjectAccessor<1, DoFHandler_type> > active_line_iterator;
    
    typedef TriaRawIterator<2,DoFCellAccessor<DoFHandler_type> >         raw_quad_iterator;
    typedef TriaIterator<2,DoFCellAccessor<DoFHandler_type> >            quad_iterator;
    typedef TriaActiveIterator<2,DoFCellAccessor<DoFHandler_type> >      active_quad_iterator;
    
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
 * @author Wolfgang Bangerth, Oliver Kayser-Herold, 1998, 2003
 */
  template <template <int> class DH>
  struct DoFIterators<DH<3> >
  {
      typedef DH<3> DoFHandler_type;
      
      typedef TriaRawIterator<3,DoFObjectAccessor<1, DoFHandler_type> >    raw_line_iterator;
      typedef TriaIterator<3,DoFObjectAccessor<1, DoFHandler_type> >       line_iterator;
      typedef TriaActiveIterator<3,DoFObjectAccessor<1, DoFHandler_type> > active_line_iterator;
      
      typedef TriaRawIterator<3,DoFObjectAccessor<2, DoFHandler_type> >    raw_quad_iterator;
      typedef TriaIterator<3,DoFObjectAccessor<2, DoFHandler_type> >       quad_iterator;
      typedef TriaActiveIterator<3,DoFObjectAccessor<2, DoFHandler_type> > active_quad_iterator;
      
      typedef TriaRawIterator<3,DoFCellAccessor<DoFHandler_type> >         raw_hex_iterator;
      typedef TriaIterator<3,DoFCellAccessor<DoFHandler_type> >            hex_iterator;
      typedef TriaActiveIterator<3,DoFCellAccessor<DoFHandler_type> >      active_hex_iterator;
      
      typedef raw_hex_iterator    raw_cell_iterator;
      typedef hex_iterator        cell_iterator;
      typedef active_hex_iterator active_cell_iterator;
      
      typedef raw_quad_iterator    raw_face_iterator;
      typedef quad_iterator        face_iterator;
      typedef active_quad_iterator active_face_iterator;    
  };

#endif // __deal2__dof_iterator_selector_h
