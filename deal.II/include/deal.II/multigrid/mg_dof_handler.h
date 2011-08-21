//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_dof_handler_h
#define __deal2__mg_dof_handler_h


#include <deal.II/base/config.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/multigrid/mg_dof_iterator_selector.h>

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
 * This class manages degrees of freedom for a multilevel hierarchy of
 * grids. It does mostly the same as does the @p DoDHandler class,
 * but it uses a separate enumeration of the degrees of freedom on
 * each level. For example, a vertex has several DoF numbers, one for
 * each level of the triangulation on which it exists.
 *
 * At present, multilevel algorithms are not fully functional, so this
 * documentation is still very brief.
 *
 * @todo This class has not yet been implemented for the use in the codimension
 * one case (<tt>spacedim != dim </tt>).
 *
//TODO:[WB] Extend MGDoFHandler doc
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim, int spacedim=dim>
class MGDoFHandler : public DoFHandler<dim,spacedim>
{
    typedef internal::MGDoFHandler::Iterators<dim,spacedim> IteratorSelector;
  public:
    typedef typename IteratorSelector::CellAccessor cell_accessor;
    typedef typename IteratorSelector::FaceAccessor face_accessor;

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
    virtual void distribute_dofs (const FiniteElement<dim,spacedim> &,
				  const unsigned int offset = 0);

				     /**
				      * Clear all data of this object
				      * and call the respective
				      * function of the base class.
				      */
    virtual void clear ();

				     /**
				      * Actually do the renumbering
				      * based on a list of new dof
				      * numbers for all the dofs.
				      *
				      * @p new_numbers is an array of
				      * integers with size equal to
				      * the number of dofs on the
				      * present level. It stores the
				      * new indices after renumbering
				      * in the order of the old
				      * indices.
				      */
    void renumber_dofs (const unsigned int               level,
			const std::vector<unsigned int> &new_numbers);

				     /**
				      * Redeclare this function of the
				      * DoFHandler basis class
				      * as otherwise it is hidden from
				      * the function with the same
				      * name, see above.
				      */
    void renumber_dofs (const std::vector<unsigned int> &new_numbers);

				     /*--------------------------------------*/

				     /**
				      *  @name Cell iterator functions
				      */
				     /*@{*/
				     /**
				      * Iterator to the first cell,
				      * used or not, on level
				      * @p level. If a level has no
				      * cells, a past-the-end iterator
				      * is returned.
				      *
				      * This function calls
				      * @p begin_raw_line in 1D and
				      * @p begin_raw_quad in 2D.
				      */
    raw_cell_iterator    begin_raw   (const unsigned int level = 0) const;

				     /**
				      * Iterator to the first used
				      * cell on level @p level.
				      *
				      * This function calls
				      * @p begin_line in 1D and
				      * @p begin_quad in 2D.
				      */
    cell_iterator        begin       (const unsigned int level = 0) const;

				     /**
				      * Iterator to the first active
				      * cell on level @p level.
				      *
				      * This function calls
				      * @p begin_active_line in 1D
				      * and @p begin_active_quad in
				      * 2D.
				      */
    active_cell_iterator begin_active(const unsigned int level = 0) const;

				     /**
				      * Iterator past the end; this
				      * iterator serves for
				      * comparisons of iterators with
				      * past-the-end or
				      * before-the-beginning states.
				      *
				      * This function calls
				      * @p end_line in 1D and
				      * @p end_quad in 2D.
				      */
    raw_cell_iterator    end () const;

				     /**
				      * Return an iterator which is
				      * the first iterator not on
				      * level. If @p level is the
				      * last level, then this returns
				      * <tt>end()</tt>.
				      */
    cell_iterator        end (const unsigned int level) const;

				     /**
				      * Return a raw iterator which is
				      * the first iterator not on
				      * level. If @p level is the
				      * last level, then this returns
				      * <tt>end()</tt>.
				      */
    raw_cell_iterator    end_raw (const unsigned int level) const;

    				     /**
				      * Return an active iterator
				      * which is the first iterator
				      * not on level. If @p level is
				      * the last level, then this
				      * returns <tt>end()</tt>.
				      */
    active_cell_iterator end_active (const unsigned int level) const;

				     /**
				      * Return an iterator pointing to
				      * the last cell, used or not.
				      *
				      * This function calls
				      * @p last_raw_line in 1D and
				      * @p last_raw_quad in 2D.
				      */
    raw_cell_iterator    last_raw () const;

				     /**
				      * Return an iterator pointing to
				      * the last cell of the level
				      * @p level, used or not.
				      *
				      * This function calls
				      * @p last_raw_line in 1D and
				      * @p last_raw_quad in 2D.
				      */
    raw_cell_iterator    last_raw (const unsigned int level) const;

				     /**
				      * Return an iterator pointing to
				      * the last used cell.
				      *
				      * This function calls
				      * @p last_line in 1D and
				      * @p last_quad in 2D.
				      */
    cell_iterator        last () const;

				     /**
				      * Return an iterator pointing to
				      * the last used cell on level
				      * @p level.
				      *
				      * This function calls
				      * @p last_line in 1D and
				      * @p last_quad in 2D.
				      */
    cell_iterator        last (const unsigned int level) const;

    				     /**
				      * Return an iterator pointing to
				      * the last active cell.
				      *
				      * This function calls
				      * @p last_active_line in 1D and
				      * @p last_active_quad in 2D.
				      */
    active_cell_iterator last_active () const;

				     /**
				      * Return an iterator pointing to
				      * the last active cell on level
				      * @p level.
				      *
				      * This function calls
				      * @p last_active_line in 1D and
				      * @p last_active_quad in 2D.
				      */
    active_cell_iterator last_active (const unsigned int level) const;
				     //@}

    				     /*---------------------------------------*/

    				     /**
				      *  @name Face iterator functions
				      */
				     /*@{*/
				     /**
				      * Iterator to the first face,
				      * used or not.
				      *
				      * This function calls
				      * @p begin_raw_line in 2D and
				      * @p begin_raw_quad in 3D.
				      */
    raw_face_iterator    begin_raw_face   () const;

				     /**
				      * Iterator to the first used
				      * face.
				      *
				      * This function calls
				      * @p begin_line in 2D and
				      * @p begin_quad in 3D.
				      */
    face_iterator        begin_face       () const;

				     /**
				      * Iterator to the first active
				      * face.
				      *
				      * This function calls
				      * @p begin_active_line in 2D
				      * and @p begin_active_quad in
				      * 3D.
				      */
    active_face_iterator begin_active_face() const;

				     /**
				      * Iterator past the end; this
				      * iterator serves for
				      * comparisons of iterators with
				      * past-the-end or
				      * before-the-beginning states.
				      *
				      * This function calls
				      * @p end_line in 2D and
				      * @p end_quad in 3D.
				      */
    raw_face_iterator    end_face () const;

				     /**
				      * Iterator past the end; this
				      * iterator serves for
				      * comparisons of iterators with
				      * past-the-end or
				      * before-the-beginning states.
				      *
				      * This is the same as
				      * <tt>end_face()</tt>.
				      */
    raw_face_iterator    end_raw_face () const;

				     /**
				      * Iterator past the end; this
				      * iterator serves for
				      * comparisons of iterators with
				      * past-the-end or
				      * before-the-beginning states.
				      *
				      * This is the same as
				      * <tt>end_face()</tt>.
				      */
    active_face_iterator end_active_face () const;

				     /**
				      * Return an iterator pointing to
				      * the last face, used or not.
				      *
				      * This function calls
				      * @p last_raw_line in 2D and
				      * @p last_raw_quad in 3D.
				      */
    raw_face_iterator    last_raw_face () const;

				     /**
				      * Return an iterator pointing to
				      * the last used face.
				      *
				      * This function calls
				      * @p last_line in 2D and
				      * @p last_quad in 3D.
				      */
    face_iterator        last_face () const;

    				     /**
				      * Return an iterator pointing to
				      * the last active face.
				      *
				      * This function calls
				      * @p last_active_line in 2D and
				      * @p last_active_quad in 3D.
				      */
    active_face_iterator last_active_face () const;
				     //@}


				     /*---------------------------------------*/

				     /**
				      * @name Line iterator functions
				      */
				     /*@{*/
				     /**
				      * Iterator to the first line,
				      * used or not, on level
				      * @p level. If a level has no
				      * lines, a past-the-end iterator
				      * is returned.
				      */
    raw_line_iterator begin_raw_line (const unsigned int level = 0) const;

				     /**
				      * Iterator to the first used
				      * line on level @p level.
				      */
    line_iterator     begin_line (const unsigned int level = 0) const;

				     /**
				      * Iterator to the first active
				      * line on level @p level.
				      */
    active_line_iterator begin_active_line(const unsigned int level = 0) const;

				     /**
				      * Iterator past the end; this
				      * iterator serves for
				      * comparisons of iterators with
				      * past-the-end or
				      * before-the-beginning states.
				      */
    raw_line_iterator end_line () const;

				     /**
				      * Return an iterator which is
				      * the first iterator not on
				      * level. If @p level is the
				      * last level, then this returns
				      * <tt>end()</tt>.
				      */
    line_iterator        end_line (const unsigned int level) const;

				     /**
				      * Return a raw iterator which is
				      * the first iterator not on
				      * level. If @p level is the
				      * last level, then this returns
				      * <tt>end()</tt>.
				      */
    raw_line_iterator    end_raw_line (const unsigned int level) const;

    				     /**
				      * Return an active iterator
				      * which is the first iterator
				      * not on level. If @p level is
				      * the last level, then this
				      * returns <tt>end()</tt>.
				      */
    active_line_iterator end_active_line (const unsigned int level) const;


				     /**
				      * Return an iterator pointing to
				      * the last line, used or not.
				      */
    raw_line_iterator    last_raw_line () const;

				     /**
				      * Return an iterator pointing to
				      * the last line of the level
				      * @p level, used or not.
				      */
    raw_line_iterator    last_raw_line (const unsigned int level) const;

				     /**
				      * Return an iterator pointing to
				      * the last used line.
				      */
    line_iterator        last_line () const;

				     /**
				      * Return an iterator pointing to
				      * the last used line on level
				      * @p level.
				      */
    line_iterator        last_line (const unsigned int level) const;

    				     /**
				      * Return an iterator pointing to
				      * the last active line.
				      */
    active_line_iterator last_active_line () const;

				     /**
				      * Return an iterator pointing to
				      * the last active line on level
				      * @p level.
				      */
    active_line_iterator last_active_line (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/

				     /**
				      *  @name Quad iterator functions*/
    				     /*@{
				      */
    				     /**
				      * Iterator to the first quad,
				      * used or not, on level
				      * @p level. If a level has no
				      * quads, a past-the-end iterator
				      * is returned.
				      */
    raw_quad_iterator    begin_raw_quad   (const unsigned int level = 0) const;

				     /**
				      * Iterator to the first used
				      * quad on level @p level.
				      */
    quad_iterator        begin_quad       (const unsigned int level = 0) const;

				     /**
				      * Iterator to the first active
				      * quad on level @p level.
				      */
    active_quad_iterator begin_active_quad(const unsigned int level = 0) const;

				     /**
				      * Iterator past the end; this
				      * iterator serves for
				      * comparisons of iterators with
				      * past-the-end or
				      * before-the-beginning states.
				      */
    raw_quad_iterator    end_quad () const;

				     /**
				      * Return an iterator which is
				      * the first iterator not on
				      * level. If @p level is the
				      * last level, then this returns
				      * <tt>end()</tt>.
				      */
    quad_iterator        end_quad (const unsigned int level) const;

				     /**
				      * Return a raw iterator which is
				      * the first iterator not on
				      * level. If @p level is the
				      * last level, then this returns
				      * <tt>end()</tt>.
				      */
    raw_quad_iterator    end_raw_quad (const unsigned int level) const;

    				     /**
				      * Return an active iterator
				      * which is the first iterator
				      * not on level. If @p level is
				      * the last level, then this
				      * returns <tt>end()</tt>.
				      */
    active_quad_iterator end_active_quad (const unsigned int level) const;


				     /**
				      * Return an iterator pointing to
				      * the last quad, used or not.
				      */
    raw_quad_iterator    last_raw_quad () const;

				     /**
				      * Return an iterator pointing to
				      * the last quad of the level
				      * @p level, used or not.
				      */
    raw_quad_iterator    last_raw_quad (const unsigned int level) const;

				     /**
				      * Return an iterator pointing to
				      * the last used quad.
				      */
    quad_iterator        last_quad () const;

				     /**
				      * Return an iterator pointing to
				      * the last used quad on level
				      * @p level.
				      */
    quad_iterator        last_quad (const unsigned int level) const;

    				     /**
				      * Return an iterator pointing to
				      * the last active quad.
				      */
    active_quad_iterator last_active_quad () const;

				     /**
				      * Return an iterator pointing to
				      * the last active quad on level
				      * @p level.
				      */
    active_quad_iterator last_active_quad (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/

				     /**
				      *  @name Hex iterator functions*/
    				     /*@{
				      */
    				     /**
				      * Iterator to the first hex,
				      * used or not, on level
				      * @p level. If a level has no
				      * hexs, a past-the-end iterator
				      * is returned.
				      */
    raw_hex_iterator    begin_raw_hex   (const unsigned int level = 0) const;

				     /**
				      * Iterator to the first used hex
				      * on level @p level.
				      */
    hex_iterator        begin_hex       (const unsigned int level = 0) const;

				     /**
				      * Iterator to the first active
				      * hex on level @p level.
				      */
    active_hex_iterator begin_active_hex(const unsigned int level = 0) const;

				     /**
				      * Iterator past the end; this
				      * iterator serves for
				      * comparisons of iterators with
				      * past-the-end or
				      * before-the-beginning states.
				      */
    raw_hex_iterator    end_hex () const;

				     /**
				      * Return an iterator which is
				      * the first iterator not on
				      * level. If @p level is the
				      * last level, then this returns
				      * <tt>end()</tt>.
				      */
    hex_iterator        end_hex (const unsigned int level) const;

				     /**
				      * Return a raw iterator which is
				      * the first iterator not on
				      * level. If @p level is the
				      * last level, then this returns
				      * <tt>end()</tt>.
				      */
    raw_hex_iterator    end_raw_hex (const unsigned int level) const;

    				     /**
				      * Return an active iterator
				      * which is the first iterator
				      * not on level. If @p level is
				      * the last level, then this
				      * returns <tt>end()</tt>.
				      */
    active_hex_iterator end_active_hex (const unsigned int level) const;


				     /**
				      * Return an iterator pointing to
				      * the last hex, used or not.
				      */
    raw_hex_iterator    last_raw_hex () const;

				     /**
				      * Return an iterator pointing to
				      * the last hex of the level
				      * @p level, used or not.
				      */
    raw_hex_iterator    last_raw_hex (const unsigned int level) const;

				     /**
				      * Return an iterator pointing to
				      * the last used hex.
				      */
    hex_iterator        last_hex () const;

				     /**
				      * Return an iterator pointing to
				      * the last used hex on level
				      * @p level.
				      */
    hex_iterator        last_hex (const unsigned int level) const;

    				     /**
				      * Return an iterator pointing to
				      * the last active hex.
				      */
    active_hex_iterator last_active_hex () const;

				     /**
				      * Return an iterator pointing to
				      * the last active hex on level
				      * @p level.
				      */
    active_hex_iterator last_active_hex (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/


    				     /**
				      * Return the number of degrees
				      * of freedom on the specified
				      * level.  Included in this
				      * number are those DoFs which
				      * are constrained by hanging
				      * nodes.
				      */
    unsigned int n_dofs (const unsigned int level) const;

				     /**
				      * Redeclare this function of the
				      * DoFHandler basis class
				      * as otherwise it is hidden from
				      * the function with the same
				      * name, see above.
				      */
    unsigned int n_dofs () const;

    				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      *
				      * This function is made virtual,
				      * since a dof handler object
				      * might be accessed through a
				      * pointers to this base class,
				      * although the actual object
				      * might be a derived class.
				      */
    virtual std::size_t memory_consumption () const;


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

  private:

				     /**
				      * We need for each vertex a list
				      * of the degree of freedom
				      * indices on each of the levels
				      * this vertex lives on. Since
				      * most vertices live only on a
				      * few levels, it is not
				      * economical to reserve indices
				      * for all the levels there are;
				      * rather, we create an object
				      * which holds the indices on
				      * those levels only where the
				      * vertex lives. To construct
				      * such an array, it is necessary
				      * to know beforehand which is
				      * the coarsest level the vertex
				      * lives on, how many levels it
				      * lives on and how many dofs
				      * there are on each vertex.  If
				      * we have this information, we
				      * can allocate exactly the
				      * amount of memory which is
				      * needed and need not handle
				      * growing arrays and the like.
				      */
    class MGVertexDoFs
    {
      public:

					 /**
					  * Constructor. This one is
					  * empty because it is
					  * difficult to make it
					  * efficient to use
					  * vector<>'s and still
					  * construct the object using
					  * the constructor. Use the
					  * @p init function to
					  * really allocate memory.
					  */
	MGVertexDoFs ();

					 /**
					  * Allocate memory and set
					  * all indices to @p -1.
					  *
					  * If @p coarsest_level is
					  * greater than @p
					  * finest_level, then no
					  * memory is allocated and
					  * the object is left in an
					  * invalid state. This is
					  * used for unused vertices.
					  */
	void init (const unsigned int coarsest_level,
		   const unsigned int finest_level,
		   const unsigned int dofs_per_vertex);

					 /**
					  * Destructor
					  */
	~MGVertexDoFs ();

					 /**
					  * Assignment operator. Will
					  * throw an exception since
					  * it can't do the work that
					  * @p init is supposed to do.
					  */
	MGVertexDoFs & operator = (const MGVertexDoFs &vertex);

					 /**
					  * Set the index with number
					  * @p dof_number of this
					  * vertex on @p level to the
					  * given index. To compute
					  * the position in the array,
					  * one has to specify how
					  * many dofs per vertex there
					  * are. It is not checked
					  * that the level number is
					  * below the number of the
					  * finest level this vertex
					  * lives on.
					  *
					  * The function is inline, so
					  * should be reasonably fast.
					  */
	void set_index (const unsigned int level,
			const unsigned int dof_number,
			const unsigned int dofs_per_vertex,
			const unsigned int index);

					 /**
					  * Return the index with
					  * number @p dof_number of
					  * this vertex on
					  * @p level. To compute the
					  * position in the array, one
					  * has to specify how many
					  * dofs per vertex there
					  * are. It is not checked
					  * that the level number is
					  * below the number of the
					  * finest level this vertex
					  * lives on.
					  *
					  * The function is inline, so
					  * should be reasonably fast.
					  */
	unsigned int get_index (const unsigned int level,
				const unsigned int dof_number,
				const unsigned int dofs_per_vertex) const;

					 /**
					  * Return the index of the
					  * coarsest level this vertex
					  * lives on.
					  */
	unsigned int get_coarsest_level () const;

					 /**
					  * Return the index of the
					  * finest level this vertex
					  * lives on.
					  */
	unsigned int get_finest_level () const;

					 /**
					  * Exception.
					  */
	DeclException0 (ExcNoMemory);
					 /**
					  * Exception.
					  */
	DeclException1 (ExcInvalidLevel,
			int,
			<< "The given level number " << arg1 << " is outside "
			<< "the range of levels this vertex lives on.");

      private:
					 /**
					  * Store the coarsest level
					  * this vertex lives on. This
					  * is used as an offset when
					  * accessing the dofs of a
					  * specific level.
					  */
	unsigned int coarsest_level;

					 /**
					  * Finest level this level
					  * lives on.  This is mostly
					  * used for error checking
					  * but can also be accessed
					  * by the function
					  * @p get_finest_level.
					  */
	unsigned int finest_level;

					 /**
					  * Array holding the indices.
					  */
	unsigned int *indices;
    };


				     /**
				      *  Return the @p i-th dof-index. This function calls
				      *  the respective function of DoFObjects.
				      */
    template <int structdim>
    unsigned int get_dof_index (const unsigned int       obj_level,
				const unsigned int       obj_index,
				const unsigned int       fe_index,
				const unsigned int       local_index) const;
				     /**
				      *  Set the @p i-th dof-index. This function calls
				      *  the respective function of DoFObjects.
				      */
    template <int structdim>
    void set_dof_index (const unsigned int       obj_level,
			const unsigned int       obj_index,
			const unsigned int       fe_index,
			const unsigned int       local_index,
			const unsigned int       global_index) const;


				     /**
				      * Reserve enough space for the
				      * MG dof indices for a given
				      * triangulation.
				      */
    void reserve_space ();

				     /**
				      * Free all used memory.
				      */
    void clear_space ();

    				     /**
				      * Space to store the DoF numbers
				      * for the different
				      * levels. Unlike the @p levels
				      * object in the
				      * DoFHandler, these are
				      * not global numbers but rather
				      * are numbers which start from
				      * zero on each level.
				      */
    std::vector<internal::DoFHandler::DoFLevel<dim>*>    mg_levels;

				     /**
				      * Space to store the DoF numbers
				      * for the faces.
				      */
    internal::DoFHandler::DoFFaces<dim> *                mg_faces;

				     /**
				      * For each vertex there is a
				      * list of indices of the degrees
				      * of freedom indices on the
				      * different levels it lives on
				      * and which are these levels.
				      */
    std::vector<MGVertexDoFs>      mg_vertex_dofs;

				     /**
				      * Vectors storing the number of
				      * degrees of freedom on each
				      * level.
				      */
    std::vector<unsigned int>      mg_used_dofs;

				     /**
				      * Make accessor objects friends.
				      */
    template <int dim1, int dim2, int dim3> friend class MGDoFAccessor;
    friend struct internal::MGDoFHandler::Implementation;
};

/*@}*/

/* ----------------------- Inline functions of MGDoFHandler -------------------*/

template <int dim, int spacedim>
inline
unsigned int MGDoFHandler<dim,spacedim>::n_dofs() const
{
  return DoFHandler<dim,spacedim>::n_dofs();
}


template <int dim, int spacedim>
inline
void MGDoFHandler<dim,spacedim>::renumber_dofs (const std::vector<unsigned int> &new_numbers)
{
  return DoFHandler<dim,spacedim>::renumber_dofs (new_numbers);
}


/* ----------------------- Inline functions of MGVertexDoFs -------------------*/

template <int dim, int spacedim>
inline
void MGDoFHandler<dim,spacedim>::MGVertexDoFs::set_index  (const unsigned int level,
						  const unsigned int dof_number,
						  const unsigned int dofs_per_vertex,
						  const unsigned int index) {
  Assert ((level >= coarsest_level) && (level <= finest_level),
	  ExcInvalidLevel(level));
  Assert (dof_number < dofs_per_vertex,
	  ExcIndexRange(dof_number, 0, dofs_per_vertex));

  indices[(level-coarsest_level)*dofs_per_vertex + dof_number] = index;
}


template <int dim, int spacedim>
inline
unsigned int
MGDoFHandler<dim,spacedim>::MGVertexDoFs::get_index  (const unsigned int level,
					     const unsigned int dof_number,
					     const unsigned int dofs_per_vertex) const {
  Assert ((level >= coarsest_level) && (level <= finest_level),
	  ExcInvalidLevel(level));
  Assert (dof_number < dofs_per_vertex,
	  ExcIndexRange (dof_number, 0, dofs_per_vertex));

  return indices[(level-coarsest_level)*dofs_per_vertex + dof_number];
}



template <>
void MGDoFHandler<1>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers);
template <>
void MGDoFHandler<2>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers);
template <>
void MGDoFHandler<3>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers);

DEAL_II_NAMESPACE_CLOSE


/*----------------------------   mg_dof.h     ---------------------------*/
#endif
/*----------------------------   mg_dof.h     ---------------------------*/
