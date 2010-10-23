//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_handler_h
#define __deal2__dof_handler_h



#include <base/config.h>
#include <base/exceptions.h>
#include <base/smartpointer.h>
#include <base/index_set.h>
#include <dofs/block_info.h>
#include <dofs/dof_iterator_selector.h>
#include <dofs/number_cache.h>
#include <dofs/function_map.h>
#include <dofs/dof_handler_policy.h>

#include <vector>
#include <map>
#include <set>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandler
  {
    template <int dim> class DoFLevel;
    template <int dim> class DoFFaces;

    struct Implementation;
  }

  namespace DoFAccessor
  {
    struct Implementation;
  }

  namespace DoFCellAccessor
  {
    struct Implementation;
  }
}


/**
 * Manage the distribution and numbering of the degrees of freedom for
 * non-multigrid algorithms.
 *
 * For each vertex, line, quad, etc, we store a list of the indices of degrees
 * of freedom living on this object. These indices refer to the unconstrained
 * degrees of freedom, i.e. constrained degrees of freedom are numbered in the
 * same way as unconstrained ones, and are only later eliminated.  This leads
 * to the fact that indices in global vectors and matrices also refer to all
 * degrees of freedom and some kind of condensation is needed to restrict the
 * systems of equations to the unconstrained degrees of freedom only. The
 * actual layout of storage of the indices is described in the internal::DoFHandler::DoFLevel class
 * documentation.
 *
 * The class offers iterators to traverse all cells, in much the same way as
 * the Triangulation class does. Using the begin() and end() functions (and
 * all their companions, like begin_active(), begin_line(), etc, just as for
 * the Triangulation class), one can obtain iterators to walk over cells, and
 * query the degree of freedom structures as well as the triangulation data.
 * These iterators are built on top of those of the Triangulation class, but
 * offer the additional information on degrees of freedom functionality than
 * pure triangulation iterators. The order in which dof iterators are
 * presented by the <tt>++</tt> and <tt>--</tt> operators is the same as that
 * for the corresponding triangulation iterators.
 *
 * The <tt>spacedim</tt> parameter has to be used if one wants to
 * solve problems in the boundary element method formulation or in an
 * equivalent one, as it is explained in the Triangulation class. If
 * not specified, this parameter takes the default value <tt>=dim</tt>
 * so that this class can be used to solve problems in the finite
 * element method formulation.
 *
 *
 * <h3>Distribution of indices for degrees of freedom</h3>
 *
 * The degrees of freedom (`dofs') are distributed on the given triangulation
 * by the function distribute_dofs(). It gets passed a finite element object
 * describing how many degrees of freedom are located on vertices, lines, etc.
 * It traverses the triangulation cell by cell and numbers the dofs of that
 * cell if not yet numbered. For non-multigrid algorithms, only active cells
 * are considered. Active cells are defined to be those cells which have no
 * children, i.e. they are the most refined ones.
 *
 * Since the triangulation is traversed starting with the cells of the coarsest
 * active level and going to more refined levels, the lowest numbers for dofs
 * are given to the largest cells as well as their bounding lines and vertices,
 * with the dofs of more refined cells getting higher numbers.
 *
 * This numbering implies very large bandwiths of the resulting matrices and
 * is thus vastly suboptimal for some solution algorithms. For this reason,
 * the DoFRenumbering class offers several algorithms to reorder the dof
 * numbering according. See there for a discussion of the implemented
 * algorithms.
 *
 *
 * <h3>Interaction with distributed meshes</h3>
 *
 * Upon construction, this class takes a reference to a triangulation
 * object. In most cases, this will be a reference to an object of
 * type Triangulation, i.e. the class that represents triangulations
 * that entirely reside on a single processor. However, it can also be
 * of type parallel::distributed::Triangulation (see, for example,
 * step-32, step-40 and in particular the @ref distributed module) in
 * which case the DoFHandler object will proceed to only manage
 * degrees of freedom on locally owned and ghost cells. This process
 * is entirely transparent to the used.
 *
 *
 * <h3>User defined renumbering schemes</h3>
 *
 * The DoFRenumbering class offers a number of renumbering schemes like the
 * Cuthill-McKey scheme. Basically, the function sets up an array in which for
 * each degree of freedom we store the new index this DoF should have after
 * renumbering. Using this array, the renumber_dofs() function of the present
 * class is called, which actually performs the change from old DoF indices to
 * the ones given in the array. In some cases, however, a user may want to
 * compute her own renumbering order; in this case, one can allocate an array
 * with one element per degree of freedom and fill it with the number that the
 * respective degree of freedom shall be assigned. This number may, for
 * example, be obtained by sorting the support points of the degrees of
 * freedom in downwind direction.  Then call the
 * <tt>renumber_dofs(vector<unsigned int>)</tt> function with the array, which
 * converts old into new degree of freedom indices.
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, int spacedim=dim>
class DoFHandler  :  public Subscriptor
{
    typedef internal::DoFHandler::Iterators<DoFHandler<dim,spacedim> > IteratorSelector;
  public:
    typedef typename IteratorSelector::CellAccessor         cell_accessor;
    typedef typename IteratorSelector::FaceAccessor         face_accessor;

    typedef typename IteratorSelector::raw_line_iterator    raw_line_iterator;
    typedef typename IteratorSelector::line_iterator        line_iterator;
    typedef typename IteratorSelector::active_line_iterator active_line_iterator;

    typedef typename IteratorSelector::raw_quad_iterator    raw_quad_iterator;
    typedef typename IteratorSelector::quad_iterator        quad_iterator;
    typedef typename IteratorSelector::active_quad_iterator active_quad_iterator;

    typedef typename IteratorSelector::raw_hex_iterator    raw_hex_iterator;
    typedef typename IteratorSelector::hex_iterator        hex_iterator;
    typedef typename IteratorSelector::active_hex_iterator active_hex_iterator;

    typedef typename IteratorSelector::raw_cell_iterator    raw_cell_iterator;
    typedef typename IteratorSelector::cell_iterator        cell_iterator;
    typedef typename IteratorSelector::active_cell_iterator active_cell_iterator;

    typedef typename IteratorSelector::raw_face_iterator    raw_face_iterator;
    typedef typename IteratorSelector::face_iterator        face_iterator;
    typedef typename IteratorSelector::active_face_iterator active_face_iterator;

				     /**
				      * Alias the @p FunctionMap type
				      * declared elsewhere.
				      */
    typedef typename dealii::FunctionMap<spacedim>::type FunctionMap;

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
				      * finite element to be used on a
				      * given cell. Since the present
				      * class only supports the same
				      * finite element to be used on
				      * all cells, the index of the
				      * finite element needs to be the
				      * same on all cells anyway, and
				      * by convention we pick zero for
				      * this value. The situation is
				      * different for hp objects
				      * (i.e. the hp::DoFHandler
				      * class) where different finite
				      * element indices may be used on
				      * different cells, and the
				      * default index there
				      * corresponds to an invalid
				      * value.
				      */
    static const unsigned int default_fe_index = 0;

				     /**
				      * Standard constructor, not
				      * initializing any data. After
				      * constructing an object with
				      * this constructor, use
				      * initialize() to make a valid
				      * DoFHandler.
				      */
    DoFHandler ();

				     /**
				      * Constructor. Take @p tria as the
				      * triangulation to work on.
				      */
    DoFHandler ( const Triangulation<dim,spacedim> &tria);

				     /**
				      * Destructor.
				      */
    virtual ~DoFHandler ();

				     /**
				      * Assign a Triangulation and a
				      * FiniteElement to the
				      * DoFHandler and compute the
				      * distribution of degrees of
				      * freedom over the mesh.
				      */
    void initialize(const Triangulation<dim,spacedim>& tria,
		    const FiniteElement<dim,spacedim>& fe);

				     /**
				      * Go through the triangulation and
				      * distribute the degrees of freedoms
				      * needed for the given finite element
				      * according to the given distribution
				      * method.
				      *
				      * The additional optional
				      * parameter @p offset allows you
				      * to reserve space for a finite
				      * number of additional vector
				      * entries in the beginning of
				      * all discretization vectors, by
				      * starting the enumeration of
				      * degrees of freedom on the grid
				      * at a nonzero value. By
				      * default, this value is of
				      * course zero.
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
    virtual void distribute_dofs (const FiniteElement<dim,spacedim> &fe,
				  const unsigned int        offset = 0);

				     /**
				      * After distribute_dofs() with
				      * an FESystem element, the block
				      * structure of global and level
				      * vectors is stored in a
				      * BlockInfo object accessible
				      * with block_info(). This
				      * function initializes the local
				      * block structure on each cell
				      * in the same object.
				      */
    void initialize_local_block_info();

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
                                        * This function is called by
                                        * the functions in
                                        * DoFRenumbering function
                                        * after computing the ordering
                                        * of the degrees of freedom.
                                        * This function is called, for
                                        * example, by the functions in
                                        * the DoFRenumbering
                                        * namespace, but it can of
                                        * course also be called from
                                        * user code.
					*
					* @arg new_number This array
					* must have a size equal to
					* the number of degrees of
					* freedom owned by the current
					* processor, i.e. the size
					* must be equal to what
					* n_locally_owned_dofs()
					* returns. If only one
					* processor participates in
					* storing the current mesh,
					* then this equals the total
					* number of degrees of
					* freedom, i.e. the result of
					* n_dofs(). The contents of
					* this array are the new
					* global indices for each
					* freedom listed in the
					* IndexSet returned by
					* locally_owned_dofs(). In the
					* case of a sequential mesh
					* this means that the array is
					* a list of new indices for
					* each of the degrees of
					* freedom on the current
					* mesh. In the case that we
					* have a
					* parallel::distributed::Triangulation
					* underlying this DoFHandler
					* object, the array is a list
					* of new indices for all the
					* locally owned degrees of
					* freedom, enumerated in the
					* same order as the currently
					* locally owned DoFs. In other
					* words, assume that degree of
					* freedom <code>i</code> is
					* currently locally owned,
					* then
					* <code>new_numbers[locally_owned_dofs().index_within_set(i)]</code>
					* returns the new global DoF
					* index of
					* <code>i</code>. Since the
					* IndexSet of
					* locally_owned_dofs() is
					* complete in the sequential
					* case, the latter convention
					* for the content of the array
					* reduces to the former in the
					* case that only one processor
					* participates in the mesh.
                                        */
    void renumber_dofs (const std::vector<unsigned int> &new_numbers);

				     /**
				      * @deprecated Use
				      * CompressedSparsityPattern instead of
				      * initializing SparsityPattern with this
				      * value, see the discussion in step-2
				      * and the @ref Sparsity module.
				      *
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
				      * matrix as well.
				      *
				      * The determination of the number of
				      * couplings can be done by simple
				      * picture drawing. An example can be
				      * found in the implementation of this
				      * function.
				      *
				      * Note that this function is most often
				      * used to determine the maximal row
				      * length for sparsity
				      * patterns. Unfortunately, while the
				      * estimates returned by this function
				      * are rather accurate in 1d and 2d, they
				      * are often significantly too high in
				      * 3d, leading the SparsityPattern class
				      * to allocate much too much memory in
				      * some cases. Unless someone comes
				      * around to improving the present
				      * function for 3d, there is not very
				      * much one can do about these cases. The
				      * typical way to work around this
				      * problem is to use an intermediate
				      * compressed sparsity pattern that only
				      * allocates memory on demand. Refer to
				      * the step-2 and step-11 example
				      * programs on how to do this. The problem
				      * is also discussed in the documentation
				      * of the module on @ref Sparsity.
				      */
    unsigned int max_couplings_between_dofs () const;

				     /**
				      * @deprecated Use
				      * CompressedSparsityPattern
				      * instead of initializing
				      * SparsityPattern with this
				      * value.
				      *
				      * Return the number of degrees of freedom
				      * located on the boundary another dof on
				      * the boundary can couple with.
				      *
				      * The number is the same as for
				      * max_couplings_between_dofs() in one
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
    raw_face_iterator    begin_raw_face   () const;

				     /**
				      *  Iterator to the first used face
				      *  on level @p level.
				      *
				      *  This function calls @p begin_line
				      *  in 2D and @p begin_quad in 3D.
				      */
    face_iterator        begin_face       () const;

				     /**
				      *  Iterator to the first active
				      *  face on level @p level.
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
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p level is
				      * the last level, then this returns
				      * <tt>end()</tt>.
				      */
    raw_face_iterator    end_raw_face () const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p level
				      * is the last level, then this returns
				      * <tt>end()</tt>.
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
				      * Access to an object informing
				      * of the block structure of the
				      * dof handler.
				      *
				      * If an FESystem is used in
				      * distribute_dofs(), degrees ofd
				      * freedom naturally split into
				      * several @ref GlossBlock
				      * "blocks". For each base element
				      * as many blocks appear as its
				      * multiplicity.
				      *
				      * At the end of
				      * distribute_dofs(), the number
				      * of degrees of freedom in each
				      * block is counted, and stored
				      * in a BlockInfo object, which
				      * can be accessed here. In an
				      * MGDoFHandler, the same is done
				      * on each level. Additionally,
				      * the block structure on each
				      * cell can be generated in this
				      * object by calling
				      * initialize_local_block_indices().
				      */
    const BlockInfo& block_info() const;


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
				      * the selected finite element
				      * object.
				      */
    const FiniteElement<dim,spacedim> & get_fe () const;

				     /**
				      * Return a constant reference to
				      * the triangulation underlying
				      * this object.
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
				      * pointers to this base class,
				      * although the actual object
				      * might be a derived class.
				      */
    virtual unsigned int memory_consumption () const;

    				     /**
				      * @todo Replace by ExcInternalError.
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
				      * Exception
				      */
    DeclException1 (ExcNewNumbersNotConsecutive,
		    int,
		    << "The given list of new dof indices is not consecutive: "
		    << "the index " << arg1 << " does not exist.");
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
				      * The object containing
				      * information on the block structure.
				      */
    BlockInfo block_info_object;

				     /**
				      * Array to store the indices for
				      * degrees of freedom located at
				      * vertices.
				      */
    std::vector<unsigned int>      vertex_dofs;



				     /**
				      * Address of the triangulation to
				      * work on.
				      */
    SmartPointer<const Triangulation<dim,spacedim>,DoFHandler<dim,spacedim> >
    tria;

				     /**
				      * Store a pointer to the finite element
				      * given latest for the distribution of
				      * dofs. In order to avoid destruction of
				      * the object before the lifetime of
				      * the DoF handler, we subscribe to
				      * the finite element object. To unlock
				      * the FE before the end of the lifetime
				      * of this DoF handler, use the <tt>clear()</tt>
				      * function (this clears all data of
				      * this object as well, though).
				      */
    SmartPointer<const FiniteElement<dim,spacedim>,DoFHandler<dim,spacedim> >
    selected_fe;

				     /**
				      * An object that describes how degrees
				      * of freedom should be distributed and
				      * renumbered.
				      */
    std_cxx1x::shared_ptr<internal::DoFHandler::Policy::PolicyBase<dim,spacedim> > policy;

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
				      * Space to store the DoF numbers
				      * for the different
				      * levels. Analogous to the
				      * <tt>levels[]</tt> tree of the
				      * Triangulation objects.
				      */
    std::vector<internal::DoFHandler::DoFLevel<dim>*> levels;

				     /**
				      * Space to store DoF numbers of
				      * faces. They are not stored in
				      * <tt>levels</tt> since faces
				      * are not organized
				      * hierarchically, but in a flat
				      * array.
				      */
    internal::DoFHandler::DoFFaces<dim> *faces;

				     /**
				      * Make accessor objects friends.
				      */
    template <int, class> friend class DoFAccessor;
    template <class> friend class DoFCellAccessor;
    friend class internal::DoFAccessor::Implementation;
    friend class internal::DoFCellAccessor::Implementation;

    friend class internal::DoFHandler::Implementation;
    friend class internal::DoFHandler::Policy::Implementation;
};




/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <> unsigned int DoFHandler<1>::n_boundary_dofs () const;
template <> unsigned int DoFHandler<1>::n_boundary_dofs (const FunctionMap &) const;
template <> unsigned int DoFHandler<1>::n_boundary_dofs (const std::set<unsigned char> &) const;


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



template <int dim, int spacedim>
inline
const FiniteElement<dim,spacedim> &
DoFHandler<dim,spacedim>::get_fe () const
{
  Assert(selected_fe!=0, ExcNotInitialized());
  return *selected_fe;
}


template <int dim, int spacedim>
inline
const Triangulation<dim,spacedim> &
DoFHandler<dim,spacedim>::get_tria () const
{
  Assert(tria != 0, ExcNotInitialized());
  return *tria;
}


template <int dim, int spacedim>
inline
const BlockInfo&
DoFHandler<dim,spacedim>::block_info () const
{
  return block_info_object;
}




#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

/*----------------------------   dof_handler.h     ---------------------------*/
#endif
/*----------------------------   dof_handler.h     ---------------------------*/
