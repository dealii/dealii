//----------------------------  dof_handler.h  ---------------------------
//    Version: $Name$
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_handler.h  ---------------------------
#ifndef __deal2__dof_handler_h
#define __deal2__dof_handler_h



#include <base/config.h>
#include <base/exceptions.h>
#include <base/smartpointer.h>
#include <dofs/function_map.h>
#include <vector>
#include <map>
#include <set>

template <int dim> class DoFCellAccessor;
template <int dim> class DoFLevel;
template<int celldim, int dim> class DoFObjectAccessor;
template <int dim> class FiniteElement;
template <int dim, typename Accessor> class TriaRawIterator;
template <int dim, typename Accessor> class TriaIterator;
template <int dim, typename Accessor> class TriaActiveIterator;
template <int dim> class Triangulation;


/**
 * Define some types which differ between the dimensions. This class
 * is analogous to the @ref{TriaDimensionInfo} class hierarchy.
 * 
 * @see DoFDimensionInfo<1>
 * @see DoFDimensionInfo<2>
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class DoFDimensionInfo
{
};



/**
 * Define some types for the DoF handling in one dimension.
 *
 * The types have the same meaning as those declared in @ref{TriaDimensionInfo<2>}.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <>
class DoFDimensionInfo<1>
{
  public:
    typedef TriaRawIterator<1,DoFCellAccessor<1> >    raw_line_iterator;
    typedef TriaIterator<1,DoFCellAccessor<1> >       line_iterator;
    typedef TriaActiveIterator<1,DoFCellAccessor<1> > active_line_iterator;

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
 * The types have the same meaning as those declared in @ref{TriaDimensionInfo<2>}.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <>
class DoFDimensionInfo<2>
{
  public:
    typedef TriaRawIterator<2,DoFObjectAccessor<1, 2> >    raw_line_iterator;
    typedef TriaIterator<2,DoFObjectAccessor<1, 2> >       line_iterator;
    typedef TriaActiveIterator<2,DoFObjectAccessor<1, 2> > active_line_iterator;
    
    typedef TriaRawIterator<2,DoFCellAccessor<2> >         raw_quad_iterator;
    typedef TriaIterator<2,DoFCellAccessor<2> >            quad_iterator;
    typedef TriaActiveIterator<2,DoFCellAccessor<2> >      active_quad_iterator;

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
 * The types have the same meaning as those declared in @ref{TriaDimensionInfo<3>}.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <>
class DoFDimensionInfo<3>
{
  public:
    typedef TriaRawIterator<3,DoFObjectAccessor<1, 3> >    raw_line_iterator;
    typedef TriaIterator<3,DoFObjectAccessor<1, 3> >       line_iterator;
    typedef TriaActiveIterator<3,DoFObjectAccessor<1, 3> > active_line_iterator;

    typedef TriaRawIterator<3,DoFObjectAccessor<2, 3> >    raw_quad_iterator;
    typedef TriaIterator<3,DoFObjectAccessor<2, 3> >       quad_iterator;
    typedef TriaActiveIterator<3,DoFObjectAccessor<2, 3> > active_quad_iterator;

    typedef TriaRawIterator<3,DoFCellAccessor<3> >         raw_hex_iterator;
    typedef TriaIterator<3,DoFCellAccessor<3> >            hex_iterator;
    typedef TriaActiveIterator<3,DoFCellAccessor<3> >      active_hex_iterator;

    typedef raw_hex_iterator    raw_cell_iterator;
    typedef hex_iterator        cell_iterator;
    typedef active_hex_iterator active_cell_iterator;

    typedef raw_quad_iterator    raw_face_iterator;
    typedef quad_iterator        face_iterator;
    typedef active_quad_iterator active_face_iterator;    
};



/**
 * Manage the distribution and numbering of the degrees of freedom for
 * non-multigrid algorithms.
 *
 * We store a list of numbers for each vertex, line, quad, etc
 * denoting the mapping between the degrees of freedom on this object
 * and the global number of this degree of freedom. The numbers refer
 * to the unconstrained degrees of freedom, i.e. constrained degrees
 * of freedom are numbered in the same way as unconstrained ones.
 * This leads to the fact that indices in global vectors and matrices
 * also refer to all degrees of freedom and some kind of condensation
 * is needed to restrict the systems of equations to the unconstrained
 * degrees of freedom only. The actual layout of storage of the indices
 * is described in the @ref{DoFLevel} class documentation.
 *
 * Finally it offers a starting point for the assemblage of the matrices
 * by offering @p{begin()} and @p{end()} functions which return iterators
 * to walk on the DoF structures as well as the triangulation data.
 * These iterators work much like those described in the documentation
 * of the @ref{Triangulation} class and of the iterator classes themselved,
 * but offer more functionality than pure triangulation iterators. The
 * order in which dof iterators are presented by the @p{++} and @p{--} operators
 * is the same as that for the alike triangulation iterators.
 *
 * This class also provides functions to create the sparsity patterns of
 * global matrices as well as matrices living on (parts of) the boundary.
 * 
 * 
 * @sect3{Distribution of indices for degrees of freedom}
 *
 * The degrees of freedom (`dofs') are distributed on the given triangulation
 * by the function @p{distribute_dofs()}. It gets passed a finite element object
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
 * the @ref{DoFRenumbering} class offers the function @p{renumber_dofs} which reorders
 * the dof numbering according to some scheme. See there for a discussion of
 * the implemented algorithms.
 *
 *
 * @sect3{User defined renumbering schemes}
 *
 * The @ref{DoFRenumbering} class offers a fixed number of renumbering
 * schemes like the Cuthill-McKey scheme. Basically, the function sets
 * up an array in which for each degree of freedom the index is stored
 * which is to be assigned by the renumbering. Using this array, the
 * @p{renumber_dofs(vector<unsigned int>)} function is called, which actually
 * does the change from old DoF indices to the ones given in the
 * array. In some cases, however, a user may want to compute her own
 * renumbering order; in this case, allocate an array with one element
 * per degree of freedom and fill it with the number that the
 * respective degree of freedom shall be assigned. This number may,
 * for example, be obtained by sorting the support points of the
 * degrees of freedom in downwind direction.  Then call the
 * @p{renumber_dofs(vector<unsigned int>)} with the array, which converts old
 * into new degree of freedom indices.
 *
 *
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class DoFHandler  :  public Subscriptor,
		     public DoFDimensionInfo<dim>
{
  public:
    typedef typename DoFDimensionInfo<dim>::raw_line_iterator raw_line_iterator;
    typedef typename DoFDimensionInfo<dim>::line_iterator line_iterator;
    typedef typename DoFDimensionInfo<dim>::active_line_iterator active_line_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_quad_iterator raw_quad_iterator;
    typedef typename DoFDimensionInfo<dim>::quad_iterator quad_iterator;
    typedef typename DoFDimensionInfo<dim>::active_quad_iterator active_quad_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_hex_iterator raw_hex_iterator;
    typedef typename DoFDimensionInfo<dim>::hex_iterator hex_iterator;
    typedef typename DoFDimensionInfo<dim>::active_hex_iterator active_hex_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_cell_iterator raw_cell_iterator;
    typedef typename DoFDimensionInfo<dim>::cell_iterator cell_iterator;
    typedef typename DoFDimensionInfo<dim>::active_cell_iterator active_cell_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_face_iterator raw_face_iterator;
    typedef typename DoFDimensionInfo<dim>::face_iterator face_iterator;
    typedef typename DoFDimensionInfo<dim>::active_face_iterator active_face_iterator;

				     /**
				      * Alias the @p{FunctionMap} type
				      * declared elsewhere.
				      */
    typedef typename FunctionMap<dim>::type FunctionMap;
    
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
    static const unsigned int invalid_dof_index = static_cast<unsigned int>(-1);
    
				     /**
				      * Constructor. Take @p{tria} as the
				      * triangulation to work on.
				      */
    DoFHandler (Triangulation<dim> &tria);
    
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
				      * The additional optional
				      * parameter @p{offset} allows you
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
				      * may want to call the @p{clear}
				      * member function which also
				      * releases the lock of this
				      * object to the finite element.
				      *
				      * This function uses the user
				      * flags of the triangulation
				      * object, so make sure you don't
				      * need them after calling this
				      * function, or if so store them.
				      */
    virtual void distribute_dofs (const FiniteElement<dim> &fe,
				  const unsigned int        offset = 0);

				     /**
				      * Clear all data of this object and
				      * especially delete the lock this object
				      * has to the finite element used the last
				      * time when @p{distribute_dofs} was called.
				      */
    virtual void clear ();
    
				     /**
				      * Actually do the renumbering based on
				      * a list of new dof numbers for all the
				      * dofs.
				      *
				      * @p{new_numbers} is an array of integers
				      * with size equal to the number of dofs
				      * on the present grid. It stores the new
				      * indices after renumbering in the
				      * order of the old indices.
				      *
				      * This function is called by the
				      * @p{renumber_dofs} function after computing
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
				      * construction of the @ref{SparsityPattern}
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
				      * @p{max_coupling_between_dofs} in one
				      * dimension less.
				      */
    unsigned int max_couplings_between_boundary_dofs () const;

				     /**
				      *  @name Cell iterator functions
				      */
				     /*@{*/
				     /**
				      *  Iterator to the first cell, used
				      *  or not, on level @p{level}. If a level
				      *  has no cells, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls @p{begin_raw_line}
				      *  in 1D and @p{begin_raw_quad} in 2D.
				      */
    raw_cell_iterator    begin_raw   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used cell
				      *  on level @p{level}.
				      *
				      *  This function calls @p{begin_line}
				      *  in 1D and @p{begin_quad} in 2D.
				      */
    cell_iterator        begin       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  cell on level @p{level}.
				      *
				      *  This function calls @p{begin_active_line}
				      *  in 1D and @p{begin_active_quad} in 2D.
				      */
    active_cell_iterator begin_active(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls @p{end_line}
				      *  in 1D and @p{end_quad} in 2D.
				      */
    raw_cell_iterator    end () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    cell_iterator        end (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_cell_iterator    end_raw (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
				      */
    active_cell_iterator end_active (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last cell, used or not.
				      *
				      *  This function calls @p{last_raw_line}
				      *  in 1D and @p{last_raw_quad} in 2D.
				      */
    raw_cell_iterator    last_raw () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  cell of the level @p{level}, used or not.
				      *
				      *  This function calls @p{last_raw_line}
				      *  in 1D and @p{last_raw_quad} in 2D.
				      */
    raw_cell_iterator    last_raw (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell.
				      *
				      *  This function calls @p{last_line}
				      *  in 1D and @p{last_quad} in 2D.
				      */
    cell_iterator        last () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell on level @p{level}.
				      *
				      *  This function calls @p{last_line}
				      *  in 1D and @p{last_quad} in 2D.
				      */
    cell_iterator        last (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active cell.
				      *
				      *  This function calls @p{last_active_line}
				      *  in 1D and @p{last_active_quad} in 2D.
				      */
    active_cell_iterator last_active () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active cell on level @p{level}.
				      *
				      *  This function calls @p{last_active_line}
				      *  in 1D and @p{last_active_quad} in 2D.
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
				      *  or not, on level @p{level}. If a level
				      *  has no faces, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls @p{begin_raw_line}
				      *  in 2D and @p{begin_raw_quad} in 3D.
				      */
    raw_face_iterator    begin_raw_face   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used face
				      *  on level @p{level}.
				      *
				      *  This function calls @p{begin_line}
				      *  in 2D and @p{begin_quad} in 3D.
				      */
    face_iterator        begin_face       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  face on level @p{level}.
				      *
				      *  This function calls @p{begin_active_line}
				      *  in 2D and @p{begin_active_quad} in 3D.
				      */
    active_face_iterator begin_active_face(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls @p{end_line}
				      *  in 2D and @p{end_quad} in 3D.
				      */
    raw_face_iterator    end_face () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    face_iterator        end_face (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_face_iterator    end_raw_face (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
				      */
    active_face_iterator end_active_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last face, used or not.
				      *
				      *  This function calls @p{last_raw_line}
				      *  in 2D and @p{last_raw_quad} in 3D.
				      */
    raw_face_iterator    last_raw_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  face of the level @p{level}, used or not.
				      *
				      *  This function calls @p{last_raw_line}
				      *  in 2D and @p{last_raw_quad} in 3D.
				      */
    raw_face_iterator    last_raw_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face.
				      *
				      *  This function calls @p{last_line}
				      *  in 2D and @p{last_quad} in 3D.
				      */
    face_iterator        last_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face on level @p{level}.
				      *
				      *  This function calls @p{last_line}
				      *  in 2D and @p{last_quad} in 3D.
				      */
    face_iterator        last_face (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active face.
				      *
				      *  This function calls @p{last_active_line}
				      *  in 2D and @p{last_active_quad} in 3D.
				      */
    active_face_iterator last_active_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active face on level @p{level}.
				      *
				      *  This function calls @p{last_active_line}
				      *  in 2D and @p{last_active_quad} in 3D.
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
				      *  or not, on level @p{level}. If a level
				      *  has no lines, a past-the-end iterator
				      *  is returned.
				      */
    raw_line_iterator    begin_raw_line   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used line
				      *  on level @p{level}.
				      */
    line_iterator        begin_line       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  line on level @p{level}.
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
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    line_iterator        end_line (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_line_iterator    end_raw_line (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
				      */
    active_line_iterator end_active_line (const unsigned int level) const;


				     /**
				      *  Return an iterator pointing to the
				      *  last line, used or not.
				      */
    raw_line_iterator    last_raw_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  line of the level @p{level}, used or not.

				      */
    raw_line_iterator    last_raw_line (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line.
				      */
    line_iterator        last_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line on level @p{level}.
				      */
    line_iterator        last_line (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active line.
				      */
    active_line_iterator last_active_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active line on level @p{level}.
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
				      *  or not, on level @p{level}. If a level
				      *  has no quads, a past-the-end iterator
				      *  is returned.
				      */
    raw_quad_iterator    begin_raw_quad   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used quad
				      *  on level @p{level}.
				      */
    quad_iterator        begin_quad       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  quad on level @p{level}.
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
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    quad_iterator        end_quad (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_quad_iterator    end_raw_quad (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
				      */
    active_quad_iterator end_active_quad (const unsigned int level) const;


				     /**
				      *  Return an iterator pointing to the
				      *  last quad, used or not.
				      */
    raw_quad_iterator    last_raw_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  quad of the level @p{level}, used or not.

				      */
    raw_quad_iterator    last_raw_quad (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad.
				      */
    quad_iterator        last_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad on level @p{level}.
				      */
    quad_iterator        last_quad (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active quad.
				      */
    active_quad_iterator last_active_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active quad on level @p{level}.
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
				      *  or not, on level @p{level}. If a level
				      *  has no hexs, a past-the-end iterator
				      *  is returned.
				      */
    raw_hex_iterator
    begin_raw_hex   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used hex
				      *  on level @p{level}.
				      */
    hex_iterator
    begin_hex       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  hex on level @p{level}.
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
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    hex_iterator        end_hex (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If @p{level} is
				      * the last level, then this returns
				      * @p{end()}.
				      */
    raw_hex_iterator    end_raw_hex (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If @p{level}
				      * is the last level, then this returns
				      * @p{end()}.
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
				      *  hex of the level @p{level}, used or not.

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
				      *  used hex on level @p{level}.
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
				      *  active hex on level @p{level}.
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
				      * reason that a @p{map} rather
				      * than a @p{set} is used is the
				      * same as descibed in the
				      * section on the
				      * @p{make_boundary_sparsity_pattern}
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
				      * Return a constant reference to the
				      * selected finite element object. This
				      * function is inline, so it should
				      * be reasonably fast.
				      */
    const FiniteElement<dim> & get_fe () const;

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

  protected:
    
				     /**
				      * Address of the triangulation to
				      * work on.
				      */
    SmartPointer<Triangulation<dim> > tria;

				     /**
				      * Store a pointer to the finite element
				      * given latest for the distribution of
				      * dofs. In order to avoid destruction of
				      * the object before the lifetime of
				      * the DoF handler, we subscribe to
				      * the finite element object. To unlock
				      * the FE before the end of the lifetime
				      * of this DoF handler, use the @p{clear()}
				      * function (this clears all data of
				      * this object as well, though).
				      */
    SmartPointer<const FiniteElement<dim> > selected_fe;

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
				      * @p{levels[]} objects to store the
				      * numbers of the degrees of freedom
				      * needed for the given element. The
				      * given element is that one which
				      * was selected when calling
				      * @p{distribute_dofs} the last time.
				      */
    void reserve_space ();

				     /**
				      * Free all used memory.
				      */
    void clear_space ();
    
				     /**
				      * Distribute dofs on the given cell,
				      * with new dofs starting with index
				      * @p{next_free_dof}. Return the next
				      * unused index number. The finite
				      * element used is the one given to
				      * @p{distribute_dofs}, which is copied
				      * to @p{selected_fe}.
				      *
				      * This function is excluded from the
				      * @p{distribute_dofs} function since
				      * it can not be implemented dimension
				      * independent.
				      */
    unsigned int distribute_dofs_on_cell (active_cell_iterator &cell,
					  unsigned int next_free_dof);

				     /**
				      * Space to store the DoF numbers for the
				      * different levels. Analogous to the
				      * @p{levels[]} tree of the @ref{Triangulation}
				      * objects.
				      */
    typename std::vector<DoFLevel<dim>*>    levels;

				     /**
				      * Store the number of dofs created last
				      * time.
				      */
    unsigned int              used_dofs;

				     /**
				      * Array to store the indices for degrees
				      * of freedom located at vertices.
				      */
    std::vector<unsigned int>      vertex_dofs;


#if (__GNUC__==2) && (__GNUC_MINOR__ < 95)
				     // this seems to be disallowed
				     // by the standard, so gcc2.95
				     // does not accept it
    friend class DoFObjectAccessor<1, dim>;
    friend class DoFObjectAccessor<2, dim>;
    friend class DoFObjectAccessor<3, dim>;
#else
				     // this, however, may grant
				     // access to too many classes,
				     // but ...
    template <int dim1, int dim2> friend class DoFObjectAccessor;
#endif
};




/* -------------- declaration of explicit specializations ------------- */

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
template <> DoFDimensionInfo<1>::raw_face_iterator DoFHandler<1>::begin_raw_face (const unsigned int) const;
template <> DoFDimensionInfo<1>::face_iterator DoFHandler<1>::begin_face (const unsigned int) const;
template <> DoFDimensionInfo<1>::active_face_iterator DoFHandler<1>::begin_active_face (const unsigned int) const;
template <> DoFDimensionInfo<1>::raw_face_iterator DoFHandler<1>::end_face () const;
template <> DoFDimensionInfo<1>::raw_face_iterator DoFHandler<1>::last_raw_face () const;
template <> DoFDimensionInfo<1>::raw_face_iterator DoFHandler<1>::last_raw_face (const unsigned int) const;
template <> DoFDimensionInfo<1>::face_iterator DoFHandler<1>::last_face () const;
template <> DoFDimensionInfo<1>::face_iterator DoFHandler<1>::last_face (const unsigned int) const;
template <> DoFDimensionInfo<1>::active_face_iterator DoFHandler<1>::last_active_face () const;
template <> DoFDimensionInfo<1>::active_face_iterator DoFHandler<1>::last_active_face (const unsigned int) const;
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
template <> DoFDimensionInfo<2>::raw_face_iterator DoFHandler<2>::begin_raw_face (const unsigned int level) const;
template <> DoFDimensionInfo<2>::face_iterator DoFHandler<2>::begin_face (const unsigned int level) const;
template <> DoFDimensionInfo<2>::active_face_iterator DoFHandler<2>::begin_active_face (const unsigned int level) const;
template <> DoFDimensionInfo<2>::raw_face_iterator DoFHandler<2>::end_face () const;
template <> DoFDimensionInfo<2>::raw_face_iterator DoFHandler<2>::last_raw_face () const;
template <> DoFDimensionInfo<2>::raw_face_iterator DoFHandler<2>::last_raw_face (const unsigned int level) const;
template <> DoFDimensionInfo<2>::face_iterator DoFHandler<2>::last_face () const;
template <> DoFDimensionInfo<2>::face_iterator DoFHandler<2>::last_face (const unsigned int level) const;
template <> DoFDimensionInfo<2>::active_face_iterator DoFHandler<2>::last_active_face () const;
template <> DoFDimensionInfo<2>::active_face_iterator DoFHandler<2>::last_active_face (const unsigned int level) const;
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

template <> unsigned int DoFHandler<1>::n_boundary_dofs () const;

template <>
unsigned int DoFHandler<1>::n_boundary_dofs (const FunctionMap &) const;
template <>
unsigned int DoFHandler<1>::n_boundary_dofs (const std::set<unsigned char> &) const;
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
unsigned int DoFHandler<dim>::n_dofs () const
{
  return used_dofs;
};



template <int dim>
inline
const FiniteElement<dim> & DoFHandler<dim>::get_fe () const
{
  Assert(selected_fe!=0, ExcNoFESelected());
  return *selected_fe;
};


/*----------------------------   dof.h     ---------------------------*/

#endif
/*----------------------------   dof.h     ---------------------------*/
