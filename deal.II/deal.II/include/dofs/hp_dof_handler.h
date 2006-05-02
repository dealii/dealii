//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------
#ifndef __deal2__hp_dof_handler_h
#define __deal2__hp_dof_handler_h



#include <base/config.h>
#include <base/exceptions.h>
#include <base/smartpointer.h>
#include <dofs/function_map.h>
#include <dofs/dof_iterator_selector.h>
#include <fe/fe_collection.h>

#include <vector>
#include <map>
#include <set>


namespace internal
{
  namespace hp
  {
    template <int> class DoFLevel;
    template <> class DoFLevel<0>;
    template <> class DoFLevel<1>;
    template <> class DoFLevel<2>;
    template <> class DoFLevel<3>;
  }
}

    

/**
 * A namespace that holds all the classes that have to do with hp finite
 * elements.
 * 
 * @ingroup hp
 */
namespace hp
{

/**
 * Manage the distribution and numbering of the degrees of freedom for
 * hp-FEM algorithms.
 *
 * @ingroup hp
 */
  template <int dim>
  class DoFHandler : public Subscriptor,
                     protected Triangulation<dim>::RefinementListener
  {
      typedef internal::DoFIteratorSelector<DoFHandler<dim> > IteratorSelector;
    public:
      typedef typename IteratorSelector::raw_line_iterator raw_line_iterator;
      typedef typename IteratorSelector::line_iterator line_iterator;
      typedef typename IteratorSelector::active_line_iterator active_line_iterator;

      typedef typename IteratorSelector::raw_quad_iterator raw_quad_iterator;
      typedef typename IteratorSelector::quad_iterator quad_iterator;
      typedef typename IteratorSelector::active_quad_iterator active_quad_iterator;

      typedef typename IteratorSelector::raw_hex_iterator raw_hex_iterator;
      typedef typename IteratorSelector::hex_iterator hex_iterator;
      typedef typename IteratorSelector::active_hex_iterator active_hex_iterator;

      typedef typename IteratorSelector::raw_cell_iterator raw_cell_iterator;
      typedef typename IteratorSelector::cell_iterator cell_iterator;
      typedef typename IteratorSelector::active_cell_iterator active_cell_iterator;

      typedef typename IteratorSelector::raw_face_iterator raw_face_iterator;
      typedef typename IteratorSelector::face_iterator face_iterator;
      typedef typename IteratorSelector::active_face_iterator active_face_iterator;

                                       /**
                                        * Alias the @p FunctionMap type
                                        * declared elsewhere.
                                        */
      typedef typename FunctionMap<dim>::type FunctionMap;
    
                                       /**
                                        * Make the dimension available
                                        * in function templates.
                                        */
      static const unsigned int dimension = dim;
    
                                       /**
                                        * When the arrays holding the
                                        * DoF indices are set up, but
                                        * before they are filled with
                                        * actual values, they are set to
                                        * an invalid value, in order to
                                        * monitor possible
                                        * problems. This invalid value
                                        * is the constant defined here.
                                        *
                                        * Please note that you should
                                        * not rely on it having a
                                        * certain value, but rather take
                                        * its symbolic name.
                                        */
      static const unsigned int invalid_dof_index = deal_II_numbers::invalid_unsigned_int;

				       /**
					* The default index of the
					* finite element to be used on
					* a given cell. For the usual,
					* non-hp ::DoFHandler class
					* that only supports the same
					* finite element to be used on
					* all cells, the index of the
					* finite element needs to be
					* the same on all cells
					* anyway, and by convention we
					* pick zero for this
					* value. The situation here is
					* different, since the hp
					* classes support the case
					* where different finite
					* element indices may be used
					* on different cells. The
					* default index consequently
					* corresponds to an invalid
					* value.
					*/
      static const unsigned int default_fe_index = deal_II_numbers::invalid_unsigned_int;

      
                                       /**
                                        * Constructor. Take @p tria as the
                                        * triangulation to work on.
                                        */
      DoFHandler (const Triangulation<dim> &tria);
    
                                       /**
                                        * Destructor.
                                        */
      virtual ~DoFHandler ();
    
                                       /**
                                        * Go through the triangulation and
                                        * distribute the degrees of freedoms
                                        * needed for the given finite element
                                        * according to the given distribution
                                        * method.
                                        *
                                        * A pointer of the transferred
                                        * finite element is
                                        * stored. Therefore, the
                                        * lifetime of the finite element
                                        * object shall be longer than
                                        * that of this object. If you
                                        * don't want this behaviour, you
                                        * may want to call the @p clear
                                        * member function which also
                                        * releases the lock of this
                                        * object to the finite element.
                                        */
      virtual void distribute_dofs (const hp::FECollection<dim> &fe);

                                       /**
                                        * Clear all data of this object and
                                        * especially delete the lock this object
                                        * has to the finite element used the last
                                        * time when @p distribute_dofs was called.
                                        */
      virtual void clear ();
    
                                       /**
                                        * Actually do the renumbering based on
                                        * a list of new dof numbers for all the
                                        * dofs.
                                        *
                                        * @p new_numbers is an array of integers
                                        * with size equal to the number of dofs
                                        * on the present grid. It stores the new
                                        * indices after renumbering in the
                                        * order of the old indices.
                                        *
                                        * This function is called by the
                                        * @p renumber_dofs function after computing
                                        * the ordering of the degrees of freedom.
                                        * However, you can call this function
                                        * yourself, which is necessary if a user
                                        * wants to implement an ordering scheme
                                        * herself, for example downwind numbering.
                                        */
      void renumber_dofs (const std::vector<unsigned int> &new_numbers);

                                       /**
                                        * Return the maximum number of
                                        * degrees of freedom a degree of freedom
                                        * in the given triangulation with the
                                        * given finite element may couple with.
                                        * This is the maximum number of entries
                                        * per line in the system matrix; this
                                        * information can therefore be used upon
                                        * construction of the SparsityPattern
                                        * object.
                                        *
                                        * The returned number is not really the
                                        * maximum number but an estimate based
                                        * on the finite element and the maximum
                                        * number of cells meeting at a vertex.
                                        * The number holds for the constrained
                                        * matrix also.
                                        *
                                        * The determination of the number of
                                        * couplings can be done by simple
                                        * picture drawing. An example can be
                                        * found in the implementation of this
                                        * function.
                                        */
      unsigned int max_couplings_between_dofs () const;

                                       /**
                                        * Return the number of degrees of freedom
                                        * located on the boundary another dof on
                                        * the boundary can couple with.
                                        *
                                        * The number is the same as for
                                        * @p max_coupling_between_dofs in one
                                        * dimension less.
                                        */
      unsigned int max_couplings_between_boundary_dofs () const;

                                       /**
                                        *  @name Cell iterator functions
                                        */
                                       /*@{*/
                                       /**
                                        *  Iterator to the first cell, used
                                        *  or not, on level @p level. If a level
                                        *  has no cells, a past-the-end iterator
                                        *  is returned.
                                        *
                                        *  This function calls @p begin_raw_line
                                        *  in 1D and @p begin_raw_quad in 2D.
                                        */
      raw_cell_iterator    begin_raw   (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first used cell
                                        *  on level @p level.
                                        *
                                        *  This function calls @p begin_line
                                        *  in 1D and @p begin_quad in 2D.
                                        */
      cell_iterator        begin       (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first active
                                        *  cell on level @p level.
                                        *
                                        *  This function calls @p begin_active_line
                                        *  in 1D and @p begin_active_quad in 2D.
                                        */
      active_cell_iterator begin_active(const unsigned int level = 0) const;

                                       /**
                                        *  Iterator past the end; this
                                        *  iterator serves for comparisons of
                                        *  iterators with past-the-end or
                                        *  before-the-beginning states.
                                        *
                                        *  This function calls @p end_line
                                        *  in 1D and @p end_quad in 2D.
                                        */
      raw_cell_iterator    end () const;

                                       /**
                                        * Return an iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      cell_iterator        end (const unsigned int level) const;
    
                                       /**
                                        * Return a raw iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      raw_cell_iterator    end_raw (const unsigned int level) const;

                                       /**
                                        * Return an active iterator which is the
                                        * first iterator not on level. If @p level
                                        * is the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      active_cell_iterator end_active (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the
                                        *  last cell, used or not.
                                        *
                                        *  This function calls @p last_raw_line
                                        *  in 1D and @p last_raw_quad in 2D.
                                        */
      raw_cell_iterator    last_raw () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  cell of the level @p level, used or not.
                                        *
                                        *  This function calls @p last_raw_line
                                        *  in 1D and @p last_raw_quad in 2D.
                                        */
      raw_cell_iterator    last_raw (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used cell.
                                        *
                                        *  This function calls @p last_line
                                        *  in 1D and @p last_quad in 2D.
                                        */
      cell_iterator        last () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used cell on level @p level.
                                        *
                                        *  This function calls @p last_line
                                        *  in 1D and @p last_quad in 2D.
                                        */
      cell_iterator        last (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active cell.
                                        *
                                        *  This function calls @p last_active_line
                                        *  in 1D and @p last_active_quad in 2D.
                                        */
      active_cell_iterator last_active () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active cell on level @p level.
                                        *
                                        *  This function calls @p last_active_line
                                        *  in 1D and @p last_active_quad in 2D.
                                        */
      active_cell_iterator last_active (const unsigned int level) const;
                                       //@}

                                       /*---------------------------------------*/

                                       /**
                                        *  @name Face iterator functions
                                        */
                                       /*@{*/
                                       /**
                                        *  Iterator to the first face, used
                                        *  or not, on level @p level. If a level
                                        *  has no faces, a past-the-end iterator
                                        *  is returned.
                                        *
                                        *  This function calls @p begin_raw_line
                                        *  in 2D and @p begin_raw_quad in 3D.
                                        */
      raw_face_iterator    begin_raw_face   (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first used face
                                        *  on level @p level.
                                        *
                                        *  This function calls @p begin_line
                                        *  in 2D and @p begin_quad in 3D.
                                        */
      face_iterator        begin_face       (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first active
                                        *  face on level @p level.
                                        *
                                        *  This function calls @p begin_active_line
                                        *  in 2D and @p begin_active_quad in 3D.
                                        */
      active_face_iterator begin_active_face(const unsigned int level = 0) const;

                                       /**
                                        *  Iterator past the end; this
                                        *  iterator serves for comparisons of
                                        *  iterators with past-the-end or
                                        *  before-the-beginning states.
                                        *
                                        *  This function calls @p end_line
                                        *  in 2D and @p end_quad in 3D.
                                        */
      raw_face_iterator    end_face () const;

                                       /**
                                        * Return an iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      face_iterator        end_face (const unsigned int level) const;
    
                                       /**
                                        * Return a raw iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      raw_face_iterator    end_raw_face (const unsigned int level) const;

                                       /**
                                        * Return an active iterator which is the
                                        * first iterator not on level. If @p level
                                        * is the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      active_face_iterator end_active_face (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the
                                        *  last face, used or not.
                                        *
                                        *  This function calls @p last_raw_line
                                        *  in 2D and @p last_raw_quad in 3D.
                                        */
      raw_face_iterator    last_raw_face () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  face of the level @p level, used or not.
                                        *
                                        *  This function calls @p last_raw_line
                                        *  in 2D and @p last_raw_quad in 3D.
                                        */
      raw_face_iterator    last_raw_face (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used face.
                                        *
                                        *  This function calls @p last_line
                                        *  in 2D and @p last_quad in 3D.
                                        */
      face_iterator        last_face () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used face on level @p level.
                                        *
                                        *  This function calls @p last_line
                                        *  in 2D and @p last_quad in 3D.
                                        */
      face_iterator        last_face (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active face.
                                        *
                                        *  This function calls @p last_active_line
                                        *  in 2D and @p last_active_quad in 3D.
                                        */
      active_face_iterator last_active_face () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active face on level @p level.
                                        *
                                        *  This function calls @p last_active_line
                                        *  in 2D and @p last_active_quad in 3D.
                                        */
      active_face_iterator last_active_face (const unsigned int level) const;
                                       //@}


                                       /*---------------------------------------*/

                                       /**
                                        *  @name Line iterator functions
                                        */
                                       /*@{*/
                                       /**
                                        *  Iterator to the first line, used
                                        *  or not, on level @p level. If a level
                                        *  has no lines, a past-the-end iterator
                                        *  is returned.
                                        */
      raw_line_iterator    begin_raw_line   (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first used line
                                        *  on level @p level.
                                        */
      line_iterator        begin_line       (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first active
                                        *  line on level @p level.
                                        */
      active_line_iterator begin_active_line(const unsigned int level = 0) const;

                                       /**
                                        *  Iterator past the end; this
                                        *  iterator serves for comparisons of
                                        *  iterators with past-the-end or
                                        *  before-the-beginning states.
                                        */
      raw_line_iterator    end_line () const;

                                       /**
                                        * Return an iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      line_iterator        end_line (const unsigned int level) const;
    
                                       /**
                                        * Return a raw iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      raw_line_iterator    end_raw_line (const unsigned int level) const;

                                       /**
                                        * Return an active iterator which is the
                                        * first iterator not on level. If @p level
                                        * is the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      active_line_iterator end_active_line (const unsigned int level) const;


                                       /**
                                        *  Return an iterator pointing to the
                                        *  last line, used or not.
                                        */
      raw_line_iterator    last_raw_line () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  line of the level @p level, used or not.

                                       */
      raw_line_iterator    last_raw_line (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used line.
                                        */
      line_iterator        last_line () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used line on level @p level.
                                        */
      line_iterator        last_line (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active line.
                                        */
      active_line_iterator last_active_line () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active line on level @p level.
                                        */
      active_line_iterator last_active_line (const unsigned int level) const;
                                       /*@}*/	  

                                       /*---------------------------------------*/

                                       /**
                                        *  @name Quad iterator functions*/
                                       /*@{
                                        */
                                       /**
                                        *  Iterator to the first quad, used
                                        *  or not, on level @p level. If a level
                                        *  has no quads, a past-the-end iterator
                                        *  is returned.
                                        */
      raw_quad_iterator    begin_raw_quad   (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first used quad
                                        *  on level @p level.
                                        */
      quad_iterator        begin_quad       (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first active
                                        *  quad on level @p level.
                                        */
      active_quad_iterator begin_active_quad(const unsigned int level = 0) const;

                                       /**
                                        *  Iterator past the end; this
                                        *  iterator serves for comparisons of
                                        *  iterators with past-the-end or
                                        *  before-the-beginning states.
                                        */
      raw_quad_iterator    end_quad () const;

                                       /**
                                        * Return an iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      quad_iterator        end_quad (const unsigned int level) const;
    
                                       /**
                                        * Return a raw iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      raw_quad_iterator    end_raw_quad (const unsigned int level) const;

                                       /**
                                        * Return an active iterator which is the
                                        * first iterator not on level. If @p level
                                        * is the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      active_quad_iterator end_active_quad (const unsigned int level) const;


                                       /**
                                        *  Return an iterator pointing to the
                                        *  last quad, used or not.
                                        */
      raw_quad_iterator    last_raw_quad () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  quad of the level @p level, used or not.

                                       */
      raw_quad_iterator    last_raw_quad (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used quad.
                                        */
      quad_iterator        last_quad () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used quad on level @p level.
                                        */
      quad_iterator        last_quad (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active quad.
                                        */
      active_quad_iterator last_active_quad () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active quad on level @p level.
                                        */
      active_quad_iterator last_active_quad (const unsigned int level) const;
                                       /*@}*/

                                       /*---------------------------------------*/

                                       /**
                                        *  @name Hex iterator functions*/
                                       /*@{
                                        */
                                       /**
                                        *  Iterator to the first hex, used
                                        *  or not, on level @p level. If a level
                                        *  has no hexs, a past-the-end iterator
                                        *  is returned.
                                        */
      raw_hex_iterator
      begin_raw_hex   (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first used hex
                                        *  on level @p level.
                                        */
      hex_iterator
      begin_hex       (const unsigned int level = 0) const;

                                       /**
                                        *  Iterator to the first active
                                        *  hex on level @p level.
                                        */
      active_hex_iterator
      begin_active_hex(const unsigned int level = 0) const;

                                       /**
                                        *  Iterator past the end; this
                                        *  iterator serves for comparisons of
                                        *  iterators with past-the-end or
                                        *  before-the-beginning states.
                                        */
      raw_hex_iterator
      end_hex () const;

                                       /**
                                        * Return an iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      hex_iterator        end_hex (const unsigned int level) const;
    
                                       /**
                                        * Return a raw iterator which is the first
                                        * iterator not on level. If @p level is
                                        * the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      raw_hex_iterator    end_raw_hex (const unsigned int level) const;

                                       /**
                                        * Return an active iterator which is the
                                        * first iterator not on level. If @p level
                                        * is the last level, then this returns
                                        * <tt>end()</tt>.
                                        */
      active_hex_iterator end_active_hex (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the
                                        *  last hex, used or not.
                                        */
      raw_hex_iterator
      last_raw_hex () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  hex of the level @p level, used or not.

                                       */
      raw_hex_iterator
      last_raw_hex (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used hex.
                                        */
      hex_iterator
      last_hex () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  used hex on level @p level.
                                        */
      hex_iterator
      last_hex (const unsigned int level) const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active hex.
                                        */
      active_hex_iterator
      last_active_hex () const;

                                       /**
                                        *  Return an iterator pointing to the last
                                        *  active hex on level @p level.
                                        */
      active_hex_iterator
      last_active_hex (const unsigned int level) const;
                                       /*@}*/

                                       /*---------------------------------------*/


                                       /**
                                        * Return number of degrees of freedom.
                                        * Included in this number are those
                                        * DoFs which are constrained by
                                        * hanging nodes.
                                        */
      unsigned int n_dofs () const;

                                       /**
                                        * Return the number of degrees of freedom
                                        * located on the boundary.
                                        */
      unsigned int n_boundary_dofs () const;

                                       /**
                                        * Return the number of degrees
                                        * of freedom located on those
                                        * parts of the boundary which
                                        * have a boundary indicator
                                        * listed in the given set. The
                                        * reason that a @p map rather
                                        * than a @p set is used is the
                                        * same as descibed in the
                                        * section on the
                                        * @p make_boundary_sparsity_pattern
                                        * function.
                                        */
      unsigned int
      n_boundary_dofs (const FunctionMap &boundary_indicators) const;

                                       /**
                                        * Same function, but with
                                        * different data type of the
                                        * argument, which is here simply
                                        * a list of the boundary
                                        * indicators under
                                        * consideration.
                                        */
      unsigned int
      n_boundary_dofs (const std::set<unsigned char> &boundary_indicators) const;

                                       /**
                                        * Return a constant reference to
                                        * the set of finite element
                                        * objects that are used by this
                                        * @p DoFHandler.
                                        */
      const hp::FECollection<dim> & get_fe () const;

                                       /**
                                        * Return a constant reference to the
                                        * triangulation underlying this object.
                                        */
      const Triangulation<dim> & get_tria () const;
    
                                       /**
                                        * Determine an estimate for the
                                        * memory consumption (in bytes)
                                        * of this object.
                                        *
                                        * This function is made virtual,
                                        * since a dof handler object
                                        * might be accessed through a
                                        * pointers to thisr base class,
                                        * although the actual object
                                        * might be a derived class.
                                        */
      virtual unsigned int memory_consumption () const;

                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcInvalidTriangulation);
                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcNoFESelected);
                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcRenumberingIncomplete);
                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcGridsDoNotMatch);
                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcInvalidBoundaryIndicator);
                                       /**
                                        * Exception
                                        */
      DeclException1 (ExcMatrixHasWrongSize,
                      int,
                      << "The matrix has the wrong dimension " << arg1);
                                       /**
                                        *  Exception
                                        */
      DeclException0 (ExcFunctionNotUseful);
                                       /**
                                        * Exception
                                        */
      DeclException1 (ExcNewNumbersNotConsecutive,
                      int,
                      << "The given list of new dof indices is not consecutive: "
                      << "the index " << arg1 << " does not exist.");
				       /**
					* Exception
					*/
      DeclException2 (ExcInvalidFEIndex,
		      int, int,
		      << "The mesh contains a cell with an active_fe_index of "
		      << arg1 << ", but the finite element collection only has "
		      << arg2 << " elements");
      
    protected:
    
                                       /**
                                        * Address of the triangulation to
                                        * work on.
                                        */
      SmartPointer<const Triangulation<dim> > tria;

                                       /**
                                        * Store a pointer to the finite
                                        * element set given latest for
                                        * the distribution of dofs. In
                                        * order to avoid destruction of
                                        * the object before the lifetime
                                        * of the DoF handler, we
                                        * subscribe to the finite
                                        * element object. To unlock the
                                        * FE before the end of the
                                        * lifetime of this DoF handler,
                                        * use the <tt>clear()</tt> function
                                        * (this clears all data of this
                                        * object as well, though).
                                        */
      SmartPointer<const hp::FECollection<dim> > finite_elements;

    private:

                                       /**
                                        * Copy constructor. I can see no reason
                                        * why someone might want to use it, so
                                        * I don't provide it. Since this class
                                        * has pointer members, making it private
                                        * prevents the compiler to provide it's
                                        * own, incorrect one if anyone chose to
                                        * copy such an object.
                                        */
      DoFHandler (const DoFHandler &);

                                       /**
                                        * Copy operator. I can see no reason
                                        * why someone might want to use it, so
                                        * I don't provide it. Since this class
                                        * has pointer members, making it private
                                        * prevents the compiler to provide it's
                                        * own, incorrect one if anyone chose to
                                        * copy such an object.
                                        */
      DoFHandler & operator = (const DoFHandler &);

                                       /**
                                        * Reserve enough space in the 
                                        * <tt>levels[]</tt> objects to store the
                                        * numbers of the degrees of freedom
                                        * needed for the given element. The
                                        * given element is that one which
                                        * was selected when calling
                                        * @p distribute_dofs the last time.
                                        */
      void reserve_space ();

				       /**
					* Do that part of reserving
					* space that pertains to
					* vertices, since this is the
					* same in all space
					* dimensions.
					*/
      void reserve_space_vertices ();
      
                                       /**
                                        * Free all used memory.
                                        */
      void clear_space ();
    
                                       /**
                                        * Distribute dofs on the given cell,
                                        * with new dofs starting with index
                                        * @p next_free_dof. Return the next
                                        * unused index number. The finite
                                        * element used is the one given to
                                        * @p distribute_dofs, which is copied
                                        * to @p selected_fe.
                                        *
                                        * This function is excluded from the
                                        * @p distribute_dofs function since
                                        * it can not be implemented dimension
                                        * independent.
                                        */
      unsigned int distribute_dofs_on_cell (active_cell_iterator &cell,
                                            unsigned int next_free_dof);


                                       /**
                                        *  Create default tables for
                                        *  the active_fe_indices in
                                        *  the
                                        *  internal::hp::DoFLevels. They
                                        *  are initialized with the a
                                        *  zero indicator, meaning
                                        *  that fe[0] is going to be
                                        *  used by default.  This
                                        *  method is called before
                                        *  refinement and before
                                        *  distribute_dofs is
                                        *  called. It ensures each
                                        *  cell has a valid
                                        *  active_fe_index.
                                        */

      void create_active_fe_table ();

                                       /**
                                        *  Methods of the
                                        *  RefinementListener coming
                                        *  from the Triangulation.
                                        *  Here they are used to
                                        *  administrate the the
                                        *  active_fe_fields during the
                                        *  spatial refinement.
                                        */
      virtual void pre_refinement_notification (const Triangulation<dim> &tria);
      virtual void post_refinement_notification (const Triangulation<dim> &tria);

				       /**
					* Set the @p local_index-th
					* degree of freedom
					* corresponding to the finite
					* element specified by @p
					* fe_index on the vertex with
					* global number @p
					* vertex_index to @p
					* global_index.
					*
					* This function is needed by
					* DoFAccessor::set_vertex_dof_index
					* when distributing degrees of
					* freedom on a mesh.
					*/
      void
      set_vertex_dof_index (const unsigned int vertex_index,
			    const unsigned int fe_index,
			    const unsigned int local_index,
			    const unsigned int global_index);

				       /**
					* Get the @p local_index-th
					* degree of freedom
					* corresponding to the finite
					* element specified by @p
					* fe_index on the vertex with
					* global number @p
					* vertex_index to @p
					* global_index.
					*
					* This function is needed by
					* DoFAccessor::vertex_dof_index,
					* which in turn is called for
					* example when doing things
					* like
					* <code>cell-@>get_dof_indices()</code>.
					*/
      unsigned int
      get_vertex_dof_index (const unsigned int vertex_index,
			    const unsigned int fe_index,
			    const unsigned int local_index) const;

                                       /**
                                        * Space to store the DoF
                                        * numbers for the different
                                        * levels. Analogous to the
                                        * <tt>levels[]</tt> tree of
                                        * the Triangulation objects.
                                        */
      std::vector<internal::hp::DoFLevel<dim>*>    levels;

                                       /**
                                        * Store the number of dofs
                                        * created last time.
                                        */
      unsigned int              used_dofs;

                                       /**
                                        * Array to store the indices
                                        * for degrees of freedom
                                        * located at vertices.
					*
					* The format used here, in the
					* form of a linked list, is
					* the same as used for the
					* arrays used in the
					* hp::DoFLevel
					* hierarchy. Starting indices
					* into this array are provided
					* by the vertex_dofs_offsets
					* field.
					*
					* Access to this field is
					* generally through the
					* get_vertex_dof_index() and
					* set_vertex_dof_index()
					* functions, encapsulating the
					* actual data format used to
					* the present class.
                                        */
      std::vector<unsigned int>      vertex_dofs;

				       /**
					* For each vertex in the
					* triangulation, store the
					* offset within the
					* vertex_dofs array where the
					* dofs for this vertex start.
					*
					* As for that array, the
					* format is the same as
					* described in the
					* documentation of
					* hp::DoFLevel.
					*
					* Access to this field is
					* generally through the
					* get_vertex_dof_index() and
					* set_vertex_dof_index()
					* functions, encapsulating the
					* actual data format used to
					* the present class.
					*/
      std::vector<unsigned int>      vertex_dofs_offsets;

                                       /**
                                        * Array to store the
                                        * information, if a cell on
                                        * some level has children or
                                        * not. It is used by the
                                        * refinement listeners as a
                                        * persistent buffer during the
                                        * refinement.
                                        */
      std::vector<std::vector<bool> *> has_children;
				      
                                       /**
                                        * Make accessor objects friends.
                                        */
      template <class DH> friend class ::DoFAccessor;

                                       /**
                                        * Make accessor objects friends.
                                        */
      template <int dim1, class DH> friend class ::DoFObjectAccessor;

                                       /**
                                        * Make accessor objects friends.
                                        */
      template <class DH> friend class ::DoFCellAccessor;

                                       /**
                                        * Likewise for DoFLevel<0>
                                        * objects since they need to
                                        * access the vertex dofs in
                                        * the functions that set and
                                        * retrieve vertex dof indices.
                                        */
      template <int> friend class internal::hp::DoFLevel;
  };




/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

  template <> unsigned int DoFHandler<1>::n_boundary_dofs () const;
  template <> unsigned int DoFHandler<1>::n_boundary_dofs (const FunctionMap &) const;
  template <> unsigned int DoFHandler<1>::n_boundary_dofs (const std::set<unsigned char> &) const;
  template <> unsigned int DoFHandler<1>::max_couplings_between_dofs () const;
  template <> unsigned int DoFHandler<1>::max_couplings_between_boundary_dofs () const;
  template <> unsigned int DoFHandler<2>::max_couplings_between_dofs () const;
  template <> unsigned int DoFHandler<2>::max_couplings_between_boundary_dofs () const;
  template <> unsigned int DoFHandler<3>::max_couplings_between_dofs () const;
  template <> unsigned int DoFHandler<3>::max_couplings_between_boundary_dofs () const;

  template <> DoFHandler<1>::raw_cell_iterator DoFHandler<1>::begin_raw (const unsigned int level) const;
  template <> DoFHandler<1>::cell_iterator DoFHandler<1>::begin (const unsigned int level) const;
  template <> DoFHandler<1>::active_cell_iterator DoFHandler<1>::begin_active (const unsigned int level) const;
  template <> DoFHandler<1>::raw_cell_iterator DoFHandler<1>::end () const;
  template <> DoFHandler<1>::raw_cell_iterator DoFHandler<1>::last_raw () const;
  template <> DoFHandler<1>::raw_cell_iterator DoFHandler<1>::last_raw (const unsigned int level) const;
  template <> DoFHandler<1>::cell_iterator DoFHandler<1>::last () const;
  template <> DoFHandler<1>::cell_iterator DoFHandler<1>::last (const unsigned int level) const;
  template <> DoFHandler<1>::active_cell_iterator DoFHandler<1>::last_active () const;
  template <> DoFHandler<1>::active_cell_iterator DoFHandler<1>::last_active (const unsigned int level) const;
  template <> DoFHandler<1>::raw_face_iterator DoFHandler<1>::begin_raw_face (const unsigned int) const;
  template <> DoFHandler<1>::face_iterator DoFHandler<1>::begin_face (const unsigned int) const;
  template <> DoFHandler<1>::active_face_iterator DoFHandler<1>::begin_active_face (const unsigned int) const;
  template <> DoFHandler<1>::raw_face_iterator DoFHandler<1>::end_face () const;
  template <> DoFHandler<1>::raw_face_iterator DoFHandler<1>::last_raw_face () const;
  template <> DoFHandler<1>::raw_face_iterator DoFHandler<1>::last_raw_face (const unsigned int) const;
  template <> DoFHandler<1>::face_iterator DoFHandler<1>::last_face () const;
  template <> DoFHandler<1>::face_iterator DoFHandler<1>::last_face (const unsigned int) const;
  template <> DoFHandler<1>::active_face_iterator DoFHandler<1>::last_active_face () const;
  template <> DoFHandler<1>::active_face_iterator DoFHandler<1>::last_active_face (const unsigned int) const;
  template <> DoFHandler<1>::raw_quad_iterator DoFHandler<1>::begin_raw_quad (const unsigned int) const;
  template <> DoFHandler<1>::quad_iterator DoFHandler<1>::begin_quad (const unsigned int) const;
  template <> DoFHandler<1>::active_quad_iterator DoFHandler<1>::begin_active_quad (const unsigned int) const;
  template <> DoFHandler<1>::raw_quad_iterator DoFHandler<1>::end_quad () const;
  template <> DoFHandler<1>::raw_quad_iterator DoFHandler<1>::last_raw_quad (const unsigned int) const;
  template <> DoFHandler<1>::quad_iterator DoFHandler<1>::last_quad (const unsigned int) const;
  template <> DoFHandler<1>::active_quad_iterator DoFHandler<1>::last_active_quad (const unsigned int) const;
  template <> DoFHandler<1>::raw_quad_iterator DoFHandler<1>::last_raw_quad () const;
  template <> DoFHandler<1>::quad_iterator DoFHandler<1>::last_quad () const;
  template <> DoFHandler<1>::active_quad_iterator DoFHandler<1>::last_active_quad () const;
  template <> DoFHandler<1>::raw_hex_iterator DoFHandler<1>::begin_raw_hex (const unsigned int) const;
  template <> DoFHandler<1>::hex_iterator DoFHandler<1>::begin_hex (const unsigned int) const;
  template <> DoFHandler<1>::active_hex_iterator DoFHandler<1>::begin_active_hex (const unsigned int) const;
  template <> DoFHandler<1>::raw_hex_iterator DoFHandler<1>::end_hex () const;
  template <> DoFHandler<1>::raw_hex_iterator DoFHandler<1>::last_raw_hex (const unsigned int) const;
  template <> DoFHandler<1>::raw_hex_iterator DoFHandler<1>::last_raw_hex () const;
  template <> DoFHandler<1>::hex_iterator DoFHandler<1>::last_hex (const unsigned int) const;
  template <> DoFHandler<1>::hex_iterator DoFHandler<1>::last_hex () const;
  template <> DoFHandler<1>::active_hex_iterator DoFHandler<1>::last_active_hex (const unsigned int) const;
  template <> DoFHandler<1>::active_hex_iterator DoFHandler<1>::last_active_hex () const;
  template <> DoFHandler<2>::raw_cell_iterator DoFHandler<2>::begin_raw (const unsigned int level) const;
  template <> DoFHandler<2>::cell_iterator DoFHandler<2>::begin (const unsigned int level) const;
  template <> DoFHandler<2>::active_cell_iterator DoFHandler<2>::begin_active (const unsigned int level) const;
  template <> DoFHandler<2>::raw_cell_iterator DoFHandler<2>::end () const;
  template <> DoFHandler<2>::raw_cell_iterator DoFHandler<2>::last_raw () const;
  template <> DoFHandler<2>::raw_cell_iterator DoFHandler<2>::last_raw (const unsigned int level) const;
  template <> DoFHandler<2>::cell_iterator DoFHandler<2>::last () const;
  template <> DoFHandler<2>::cell_iterator DoFHandler<2>::last (const unsigned int level) const;
  template <> DoFHandler<2>::active_cell_iterator DoFHandler<2>::last_active () const;
  template <> DoFHandler<2>::active_cell_iterator DoFHandler<2>::last_active (const unsigned int level) const;
  template <> DoFHandler<2>::raw_face_iterator DoFHandler<2>::begin_raw_face (const unsigned int level) const;
  template <> DoFHandler<2>::face_iterator DoFHandler<2>::begin_face (const unsigned int level) const;
  template <> DoFHandler<2>::active_face_iterator DoFHandler<2>::begin_active_face (const unsigned int level) const;
  template <> DoFHandler<2>::raw_face_iterator DoFHandler<2>::end_face () const;
  template <> DoFHandler<2>::raw_face_iterator DoFHandler<2>::last_raw_face () const;
  template <> DoFHandler<2>::raw_face_iterator DoFHandler<2>::last_raw_face (const unsigned int level) const;
  template <> DoFHandler<2>::face_iterator DoFHandler<2>::last_face () const;
  template <> DoFHandler<2>::face_iterator DoFHandler<2>::last_face (const unsigned int level) const;
  template <> DoFHandler<2>::active_face_iterator DoFHandler<2>::last_active_face () const;
  template <> DoFHandler<2>::active_face_iterator DoFHandler<2>::last_active_face (const unsigned int level) const;
  template <> DoFHandler<2>::raw_hex_iterator DoFHandler<2>::begin_raw_hex (const unsigned int) const;
  template <> DoFHandler<2>::hex_iterator DoFHandler<2>::begin_hex (const unsigned int) const;
  template <> DoFHandler<2>::active_hex_iterator DoFHandler<2>::begin_active_hex (const unsigned int) const;
  template <> DoFHandler<2>::raw_hex_iterator DoFHandler<2>::end_hex () const;
  template <> DoFHandler<2>::raw_hex_iterator DoFHandler<2>::last_raw_hex (const unsigned int) const;
  template <> DoFHandler<2>::raw_hex_iterator DoFHandler<2>::last_raw_hex () const;
  template <> DoFHandler<2>::hex_iterator DoFHandler<2>::last_hex (const unsigned int) const;
  template <> DoFHandler<2>::hex_iterator DoFHandler<2>::last_hex () const;
  template <> DoFHandler<2>::active_hex_iterator DoFHandler<2>::last_active_hex (const unsigned int) const;
  template <> DoFHandler<2>::active_hex_iterator DoFHandler<2>::last_active_hex () const;
  template <> DoFHandler<3>::raw_cell_iterator DoFHandler<3>::begin_raw (const unsigned int level) const;
  template <> DoFHandler<3>::cell_iterator DoFHandler<3>::begin (const unsigned int level) const;
  template <> DoFHandler<3>::active_cell_iterator DoFHandler<3>::begin_active (const unsigned int level) const;
  template <> DoFHandler<3>::raw_cell_iterator DoFHandler<3>::end () const;
  template <> DoFHandler<3>::raw_cell_iterator DoFHandler<3>::last_raw () const;
  template <> DoFHandler<3>::raw_cell_iterator DoFHandler<3>::last_raw (const unsigned int level) const;
  template <> DoFHandler<3>::cell_iterator DoFHandler<3>::last () const;
  template <> DoFHandler<3>::cell_iterator DoFHandler<3>::last (const unsigned int level) const;
  template <> DoFHandler<3>::active_cell_iterator DoFHandler<3>::last_active () const;
  template <> DoFHandler<3>::active_cell_iterator DoFHandler<3>::last_active (const unsigned int level) const;
  template <> DoFHandler<3>::raw_face_iterator DoFHandler<3>::begin_raw_face (const unsigned int level) const;
  template <> DoFHandler<3>::face_iterator DoFHandler<3>::begin_face (const unsigned int level) const;
  template <> DoFHandler<3>::active_face_iterator DoFHandler<3>::begin_active_face (const unsigned int level) const;
  template <> DoFHandler<3>::raw_face_iterator DoFHandler<3>::end_face () const;
  template <> DoFHandler<3>::raw_face_iterator DoFHandler<3>::last_raw_face () const;
  template <> DoFHandler<3>::raw_face_iterator DoFHandler<3>::last_raw_face (const unsigned int level) const;
  template <> DoFHandler<3>::face_iterator DoFHandler<3>::last_face () const;
  template <> DoFHandler<3>::face_iterator DoFHandler<3>::last_face (const unsigned int level) const;
  template <> DoFHandler<3>::active_face_iterator DoFHandler<3>::last_active_face () const;
  template <> DoFHandler<3>::active_face_iterator DoFHandler<3>::last_active_face (const unsigned int level) const;

  template <>
  unsigned int DoFHandler<1>::distribute_dofs_on_cell (active_cell_iterator &cell,
                                                         unsigned int          next_free_dof);
  template <>
  unsigned int DoFHandler<2>::distribute_dofs_on_cell (active_cell_iterator &cell,
                                                         unsigned int          next_free_dof);
  template <>
  unsigned int DoFHandler<3>::distribute_dofs_on_cell (active_cell_iterator &cell,
                                                         unsigned int          next_free_dof);
  template <> void DoFHandler<1>::renumber_dofs (const std::vector<unsigned int> &new_numbers);
  template <> void DoFHandler<2>::renumber_dofs (const std::vector<unsigned int> &new_numbers);
  template <> void DoFHandler<3>::renumber_dofs (const std::vector<unsigned int> &new_numbers);
  template <> void DoFHandler<1>::reserve_space ();
  template <> void DoFHandler<2>::reserve_space ();
  template <> void DoFHandler<3>::reserve_space ();

/* ----------------------- Inline functions ---------------------------------- */

  template <int dim>
  inline
  unsigned int
  DoFHandler<dim>::n_dofs () const
  {
    return used_dofs;
  }



  template <int dim>
  inline
  const hp::FECollection<dim> &
  DoFHandler<dim>::get_fe () const
  {
    return *finite_elements;
  }


  template <int dim>
  inline
  const Triangulation<dim> &
  DoFHandler<dim>::get_tria () const
  {
    return *tria;
  }



  template <int dim>
  inline
  unsigned int
  DoFHandler<dim>::
  get_vertex_dof_index (const unsigned int vertex_index,
			const unsigned int fe_index,
			const unsigned int local_index) const
  {
    Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
	    ExcMessage ("You need to specify a FE index when working "
			"with hp DoFHandlers"));
    Assert (finite_elements != 0,
	    ExcMessage ("No finite element collection is associated with "
			"this DoFHandler"));
    Assert (local_index < (*finite_elements)[fe_index].dofs_per_vertex,
	    ExcIndexRange(local_index, 0,
			  (*finite_elements)[fe_index].dofs_per_vertex));

				     // hop along the list of index
				     // sets until we find the one
				     // with the correct fe_index, and
				     // then poke into that
				     // part. trigger an exception if
				     // we can't find a set for this
				     // particular fe_index
    const unsigned int starting_offset = vertex_dofs_offsets[vertex_index];
    const unsigned int *pointer        = &vertex_dofs[starting_offset];
    while (true)
      {
	Assert (*pointer != deal_II_numbers::invalid_unsigned_int,
		ExcInternalError());
	if (*pointer == fe_index)
	  return *(pointer + 1 + local_index);
	else
	  pointer += (*finite_elements)[*pointer].dofs_per_vertex;
      }
  }  



  template <int dim>
  inline
  void
  DoFHandler<dim>::
  set_vertex_dof_index (const unsigned int vertex_index,
			const unsigned int fe_index,
			const unsigned int local_index,
			const unsigned int global_index)
  {
    Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
	    ExcMessage ("You need to specify a FE index when working "
			"with hp DoFHandlers"));
    Assert (finite_elements != 0,
	    ExcMessage ("No finite element collection is associated with "
			"this DoFHandler"));
    Assert (local_index < (*finite_elements)[fe_index].dofs_per_vertex,
	    ExcIndexRange(local_index, 0,
			  (*finite_elements)[fe_index].dofs_per_vertex));

				     // hop along the list of index
				     // sets until we find the one
				     // with the correct fe_index, and
				     // then poke into that
				     // part. trigger an exception if
				     // we can't find a set for this
				     // particular fe_index
    const unsigned int starting_offset = vertex_dofs_offsets[vertex_index];
    unsigned int *pointer              = &vertex_dofs[starting_offset];
    while (true)
      {
	Assert (*pointer != deal_II_numbers::invalid_unsigned_int,
		ExcInternalError());
	if (*pointer == fe_index)
	  {
	    *(pointer + 1 + local_index) = global_index;
	    return;
	  }
	else
	  pointer += (*finite_elements)[*pointer].dofs_per_vertex;
      }
  }  
  
#endif
    
}

#endif

  
