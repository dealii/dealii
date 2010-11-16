//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
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
#include <base/template_constraints.h>
#include <base/smartpointer.h>
#include <dofs/function_map.h>
#include <dofs/dof_iterator_selector.h>
#include <dofs/number_cache.h>
#include <hp/fe_collection.h>

#include <vector>
#include <map>
#include <set>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace hp
  {
    template <int> class DoFLevel;
    template <int> class DoFFaces;
    template <int> class DoFObjects;

    namespace DoFHandler
    {
      struct Implementation;
    }
  }
}

namespace internal
{
  namespace DoFAccessor
  {
    struct Implementation;
  }

  namespace DoFCellAccessor
  {
    struct Implementation;
  }
}



namespace hp
{

/**
 * Manage the distribution and numbering of the degrees of freedom for
 * hp-FEM algorithms.
 *
 * This class has not yet been implemented for the use in the codimension
 * one case (<tt>spacedim != dim </tt>).
 *
 * @ingroup dofs
 * @ingroup hp
 */
  template <int dim, int spacedim=dim>
  class DoFHandler : public Subscriptor,
                     protected Triangulation<dim,spacedim>::RefinementListener
  {
      typedef internal::DoFHandler::Iterators<DoFHandler<dim,spacedim> > IteratorSelector;
    public:
      typedef typename IteratorSelector::CellAccessor         cell_accessor;
      typedef typename IteratorSelector::FaceAccessor         face_accessor;

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
      typedef typename FunctionMap<spacedim>::type FunctionMap;

                                       /**
                                        * Make the dimension available
                                        * in function templates.
                                        */
      static const unsigned int dimension = dim;

                                       /**
                                        * Make the space dimension available
                                        * in function templates.
                                        */
      static const unsigned int space_dimension = spacedim;

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
      static const unsigned int invalid_dof_index = numbers::invalid_unsigned_int;

				       /**
					* The default index of the
					* finite element to be used on
					* a given cell. For the usual,
					* non-hp dealii::DoFHandler class
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
      static const unsigned int default_fe_index = numbers::invalid_unsigned_int;


                                       /**
                                        * Constructor. Take @p tria as the
                                        * triangulation to work on.
                                        */
      DoFHandler (const Triangulation<dim,spacedim> &tria);

                                       /**
                                        * Destructor.
                                        */
      virtual ~DoFHandler ();

                                       /**
                                        * Go through the triangulation and
                                        * distribute the degrees of freedoms
                                        * needed for the given finite element
                                        * according to the current distribution
                                        * of active fe indices.
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
      virtual void distribute_dofs (const hp::FECollection<dim,spacedim> &fe);

                                       /**
                                        * Go through the triangulation and set
                                        * the active FE indices of all active
                                        * cells to the values given in @p
                                        * active_fe_indices.
                                        */
      void set_active_fe_indices (const std::vector<unsigned int>& active_fe_indices);

                                       /**
                                        * Go through the triangulation and
                                        * store the active FE indices of all
                                        * active cells to the vector @p
                                        * active_fe_indices. This vector is
                                        * resized, if necessary.
                                        */
      void get_active_fe_indices (std::vector<unsigned int>& active_fe_indices) const;

                                       /**
                                        * Clear all data of this object and
                                        * especially delete the lock this object
                                        * has to the finite element used the last
                                        * time when @p distribute_dofs was called.
                                        */
      virtual void clear ();

                                       /**
                                        * Renumber degrees of freedom based on
                                        * a list of new dof numbers for all the
                                        * dofs.
                                        *
                                        * @p new_numbers is an array of integers
                                        * with size equal to the number of dofs
                                        * on the present grid. It stores the new
                                        * indices after renumbering in the
                                        * order of the old indices.
                                        *
                                        * This function is called by
                                        * the functions in
                                        * DoFRenumbering function
                                        * after computing the ordering
                                        * of the degrees of freedom.
                                        * However, you can call this
                                        * function yourself, which is
                                        * necessary if a user wants to
                                        * implement an ordering scheme
                                        * herself, for example
                                        * downwind numbering.
					*
					* The @p new_number array must
					* have a size equal to the
					* number of degrees of
					* freedom. Each entry must
					* state the new global DoF
					* number of the degree of
					* freedom referenced.
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
                                        * As for
                                        * DoFHandler::max_couplings_between_dofs(),
                                        * the result of this function is often
                                        * not very accurate for 3d and/or high
                                        * polynomial degrees. The consequences
                                        * are discussed in the documentation
                                        * of the module on @ref Sparsity.
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
					*  or not.
					*
					*  This function calls @p begin_raw_line
					*  in 2D and @p begin_raw_quad in 3D.
					*/
      raw_face_iterator    begin_raw_face   () const;

				       /**
					*  Iterator to the first used face.
					*
					*  This function calls @p begin_line
					*  in 2D and @p begin_quad in 3D.
					*/
      face_iterator        begin_face       () const;

				       /**
					*  Iterator to the first active face.
					*
					*  This function calls @p begin_active_line
					*  in 2D and @p begin_active_quad in 3D.
					*/
      active_face_iterator begin_active_face() const;

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
					* Return a raw iterator past-the-end.
					* This is the same as <tt>end_face()</tt>.
					*/
      raw_face_iterator    end_raw_face () const;

				       /**
					* Return an active iterator past-the-end.
					* This is the same as <tt>end_face()</tt>.
					*/
      active_face_iterator end_active_face () const;

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
					*  used face.
					*
					*  This function calls @p last_line
					*  in 2D and @p last_quad in 3D.
					*/
      face_iterator        last_face () const;

				       /**
					*  Return an iterator pointing to the last
					*  active face.
					*
					*  This function calls @p last_active_line
					*  in 2D and @p last_active_quad in 3D.
					*/
      active_face_iterator last_active_face () const;

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
					* Return the global number of
					* degrees of freedom. If the
					* current object handles all
					* degrees of freedom itself
					* (even if you may intend to
					* solve your linear system in
					* parallel, such as in step-17
					* or step-18), then this number
					* equals the number of locally
					* owned degrees of freedom since
					* this object doesn't know
					* anything about what you want
					* to do with it and believes
					* that it owns every degree of
					* freedom it knows about.
					*
					* On the other hand, if this
					* object operates on a
					* parallel::distributed::Triangulation
					* object, then this function
					* returns the global number of
					* degrees of freedom,
					* accumulated over all
					* processors.
					*
					* In either case, included in
					* the returned number are those
					* DoFs which are constrained by
					* hanging nodes, see @ref constraints.
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
					* Return the number of
					* degrees of freedom that
					* belong to this
					* process.
					*
					* If this is a sequential job,
					* then the result equals that
					* produced by n_dofs(). On the
					* other hand, if we are
					* operating on a
					* parallel::distributed::Triangulation,
					* then it includes only the
					* degrees of freedom that the
					* current processor owns. Note
					* that in this case this does
					* not include all degrees of
					* freedom that have been
					* distributed on the current
					* processor's image of the mesh:
					* in particular, some of the
					* degrees of freedom on the
					* interface between the cells
					* owned by this processor and
					* cells owned by other
					* processors may be theirs, and
					* degrees of freedom on ghost
					* cells are also not necessarily
					* included.
					*/
      unsigned int n_locally_owned_dofs() const;

				       /**
					* Return an IndexSet describing
					* the set of locally owned DoFs
					* as a subset of
					* 0..n_dofs(). The number of
					* elements of this set equals
					* n_locally_owned_dofs().
					*/
      const IndexSet & locally_owned_dofs() const;


				       /**
					* Returns a vector that
					* stores the locally owned
					* DoFs of each processor. If
					* you are only interested in
					* the number of elements
					* each processor owns then
					* n_dofs_per_processor() is
					* a better choice.
					*
					* If this is a sequential job,
					* then the vector has a single
					* element that equals the
					* IndexSet representing the
					* entire range [0,n_dofs()].
					*/
      const std::vector<IndexSet> &
      locally_owned_dofs_per_processor () const;

				       /**
					* Return a vector that
					* stores the number of
					* degrees of freedom each
					* processor that
					* participates in this
					* triangulation owns
					* locally. The sum of all
					* these numbers equals the
					* number of degrees of
					* freedom that exist
					* globally, i.e. what
					* n_dofs() returns.
					*
					* Each element of the vector
					* returned by this function
					* equals the number of
					* elements of the
					* corresponding sets
					* returned by
					* global_dof_indices().
					*
					* If this is a sequential job,
					* then the vector has a single
					* element equal to n_dofs().
					*/
      const std::vector<unsigned int> &
      n_locally_owned_dofs_per_processor () const;

                                       /**
                                        * Return a constant reference to
                                        * the set of finite element
                                        * objects that are used by this
                                        * @p DoFHandler.
                                        */
      const hp::FECollection<dim,spacedim> & get_fe () const;

                                       /**
                                        * Return a constant reference to the
                                        * triangulation underlying this object.
                                        */
      const Triangulation<dim,spacedim> & get_tria () const;

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
				       /**
					*  Exception
					*/
      DeclException1 (ExcInvalidLevel,
		      int,
		      << "The given level " << arg1
		      << " is not in the valid range!");
				       /**
					* Exception
					*/
      DeclException0 (ExcFacesHaveNoLevel);
				       /**
					* The triangulation level you
					* accessed is empty.
					*/
      DeclException1 (ExcEmptyLevel,
		      int,
		      << "You tried to do something on level " << arg1
		      << ", but this level is empty.");

    protected:

                                       /**
                                        * Address of the triangulation to
                                        * work on.
                                        */
      SmartPointer<const Triangulation<dim,spacedim>,DoFHandler<dim,spacedim> > tria;

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
      SmartPointer<const hp::FECollection<dim,spacedim>,hp::DoFHandler<dim,spacedim> > finite_elements;

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
                                        * Free all used memory.
                                        */
      void clear_space ();

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
      virtual void pre_refinement_notification (const Triangulation<dim,spacedim> &tria);
      virtual void post_refinement_notification (const Triangulation<dim,spacedim> &tria);

				       /**
					* Compute identities between
					* DoFs located on
					* vertices. Called from
					* distribute_dofs().
					*/
      void
      compute_vertex_dof_identities (std::vector<unsigned int> &new_dof_indices) const;

				       /**
					* Compute identities between
					* DoFs located on
					* lines. Called from
					* distribute_dofs().
					*/
      void
      compute_line_dof_identities (std::vector<unsigned int> &new_dof_indices) const;

				       /**
					* Compute identities between
					* DoFs located on
					* quads. Called from
					* distribute_dofs().
					*/
      void
      compute_quad_dof_identities (std::vector<unsigned int> &new_dof_indices) const;

				       /**
					* Renumber the objects with
					* the given and all lower
					* structural dimensions,
					* i.e. renumber vertices by
					* giving a template argument
					* of zero to the int2type
					* argument, lines and vertices
					* with one, etc.
					*
					* Note that in contrast to the
					* public renumber_dofs()
					* function, these internal
					* functions do not ensure that
					* the new DoFs are
					* contiguously numbered. The
					* function may therefore also
					* be used to assign different
					* DoFs the same number, for
					* example to unify hp DoFs
					* corresponding to different
					* finite elements but
					* co-located on the same
					* entity.
					*/
      void renumber_dofs_internal (const std::vector<unsigned int> &new_numbers,
				   internal::int2type<0>);

      void renumber_dofs_internal (const std::vector<unsigned int> &new_numbers,
				   internal::int2type<1>);

      void renumber_dofs_internal (const std::vector<unsigned int> &new_numbers,
				   internal::int2type<2>);

      void renumber_dofs_internal (const std::vector<unsigned int> &new_numbers,
				   internal::int2type<3>);

                                       /**
                                        * Space to store the DoF
                                        * numbers for the different
                                        * levels. Analogous to the
                                        * <tt>levels[]</tt> tree of
                                        * the Triangulation objects.
                                        */
      std::vector<internal::hp::DoFLevel<dim>*>    levels;
                                       /**
                                        * Space to store the DoF
                                        * numbers for the faces.
                                        * Analogous to the
                                        * <tt>faces</tt> pointer of
                                        * the Triangulation objects.
                                        */
      internal::hp::DoFFaces<dim> * faces;

				       /**
					* A structure that contains all
					* sorts of numbers that
					* characterize the degrees of
					* freedom this object works on.
					*
					* For most members of this
					* structure, there is an
					* accessor function in this
					* class that returns its value.
					*/
      internal::DoFHandler::NumberCache number_cache;

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
					* DoFAccessor::get_vertex_dof_index() and
					* DoFAccessor::set_vertex_dof_index()
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
					* Accessor::get_vertex_dof_index() and
					* Accessor::set_vertex_dof_index()
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
      template <int, class> friend class dealii::DoFAccessor;
      template <class> friend class dealii::DoFCellAccessor;
      friend class internal::DoFAccessor::Implementation;
      friend class internal::DoFCellAccessor::Implementation;

                                       /**
                                        * Likewise for DoFLevel
                                        * objects since they need to
                                        * access the vertex dofs in
                                        * the functions that set and
                                        * retrieve vertex dof indices.
                                        */
      template <int> friend class internal::hp::DoFLevel;
      template <int> friend class internal::hp::DoFObjects;
      friend class internal::hp::DoFHandler::Implementation;
  };



#ifndef DOXYGEN


/* ----------------------- Inline functions ---------------------------------- */

  template <int dim, int spacedim>
  inline
  unsigned int
  DoFHandler<dim,spacedim>::n_dofs () const
  {
    return number_cache.n_global_dofs;
  }


  template <int dim, int spacedim>
  unsigned int
  DoFHandler<dim, spacedim>::n_locally_owned_dofs() const
  {
    return number_cache.n_locally_owned_dofs;
  }


  template <int dim, int spacedim>
  const IndexSet &
  DoFHandler<dim, spacedim>::locally_owned_dofs() const
  {
    return number_cache.locally_owned_dofs;
  }


  template <int dim, int spacedim>
  const std::vector<unsigned int> &
  DoFHandler<dim, spacedim>::n_locally_owned_dofs_per_processor() const
  {
    return number_cache.n_locally_owned_dofs_per_processor;
  }


  template <int dim, int spacedim>
  const std::vector<IndexSet> &
  DoFHandler<dim, spacedim>::locally_owned_dofs_per_processor () const
  {
    return number_cache.locally_owned_dofs_per_processor;
  }



  template<int dim, int spacedim>
  inline
  const hp::FECollection<dim,spacedim> &
  DoFHandler<dim,spacedim>::get_fe () const
  {
    Assert (finite_elements != 0,
	    ExcMessage ("No finite element collection is associated with "
			"this DoFHandler"));
    return *finite_elements;
  }


  template<int dim, int spacedim>
  inline
  const Triangulation<dim,spacedim> &
  DoFHandler<dim,spacedim>::get_tria () const
  {
    return *tria;
  }




#endif

}

DEAL_II_NAMESPACE_CLOSE

#endif


