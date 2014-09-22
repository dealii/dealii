// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#ifndef __deal2__mg_dof_handler_h
#define __deal2__mg_dof_handler_h

//TODO: MgDoFHandler should be marked deprecated, but is used everywhere...
//#warning MGDoFHandler is deprecated

#include <deal.II/base/config.h>
#include <deal.II/dofs/dof_handler.h>
//#include <deal.II/multigrid/mg_dof_iterator_selector.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MGDoFHandler
  {
    struct Implementation;
  }
}


/*!@addtogroup mg */
/*@{*/


/**
 * @deprecated All functionality of this class has been moved to
 * DoFHandler. Thus, this class is only a wrapper left in the library
 * for compatibility reasons.
 *
 * This class manages degrees of freedom for a multilevel hierarchy of
 * grids. It does mostly the same as does the @p DoDHandler class,
 * but it uses a separate enumeration of the degrees of freedom on
 * each level. For example, a vertex has several DoF numbers, one for
 * each level of the triangulation on which it exists.
 *
 * @todo This class has not yet been implemented for the use in the codimension
 * one case (<tt>spacedim != dim </tt>).
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, 1998, 1999 Markus Bürg, Timo Heister, Guido Kanschat, 2012
 */
template <int dim, int spacedim=dim>
class MGDoFHandler : public DoFHandler<dim,spacedim>
{
  typedef dealii::internal::DoFHandler::Iterators<DoFHandler<dim,spacedim>, true> IteratorSelector;
public:
  typedef typename IteratorSelector::CellAccessor cell_accessor;
  typedef typename IteratorSelector::FaceAccessor face_accessor;

  typedef typename IteratorSelector::raw_line_iterator raw_line_iterator;
  typedef typename IteratorSelector::line_iterator line_iterator;
  typedef typename IteratorSelector::active_line_iterator active_line_iterator;

  typedef typename IteratorSelector::quad_iterator quad_iterator;
  typedef typename IteratorSelector::active_quad_iterator active_quad_iterator;

  typedef typename IteratorSelector::hex_iterator hex_iterator;
  typedef typename IteratorSelector::active_hex_iterator active_hex_iterator;

  typedef typename IteratorSelector::cell_iterator cell_iterator;
  typedef typename IteratorSelector::active_cell_iterator active_cell_iterator;

  typedef typename IteratorSelector::face_iterator face_iterator;
  typedef typename IteratorSelector::active_face_iterator active_face_iterator;

  /**
   * Make the dimension and the space_dimension available
   * in function templates.
   */
  static const unsigned int dimension = dim;

  static const unsigned int space_dimension = spacedim;

  /**
   * Default constructor, which
   * will require a call to
   * initialize() later to set the Triangulation.
   */
  MGDoFHandler ();

  /**
   * Constructor. Take @p tria as
   * the triangulation to work on.
   */
  MGDoFHandler (const Triangulation<dim, spacedim> &tria);

  /**
   * Destructor
   */
  virtual ~MGDoFHandler ();

  /**
   * Go through the triangulation
   * and distribute the degrees of
   * freedoms needed for the given
   * finite element according to
   * the given distribution
   * method. We first call the
   * DoFHandler's function
   * and then distribute the
   * levelwise numbers.
   *
   * A copy of the transferred
   * finite element is stored.
   */
  virtual void distribute_dofs (const FiniteElement<dim,spacedim> &);

  /**
   *  @name Cell iterator functions
   */
  /*@{*/
  /**
   * Iterator to the first used
   * cell on level @p level.
   */
  cell_iterator        begin       (const unsigned int level = 0) const;

  /**
   * Iterator past the end; this
   * iterator serves for
   * comparisons of iterators with
   * past-the-end or
   * before-the-beginning states.
   */
  cell_iterator        end () const;

  /**
   * Return an iterator which is
   * the first iterator not on
   * level. If @p level is the
   * last level, then this returns
   * <tt>end()</tt>.
   */
  cell_iterator        end (const unsigned int level) const;

  //@}

} DEAL_II_DEPRECATED;

/*@}*/


DEAL_II_NAMESPACE_CLOSE


#endif
