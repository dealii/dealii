//----------------------------  mg_dof_handler.h  ---------------------------
//    Version: $Name$
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_dof_handler.h  ---------------------------
#ifndef __deal2__mg_dof_handler_h
#define __deal2__mg_dof_handler_h


/*----------------------------   mg_dof.h     ---------------------------*/


#include <dofs/dof_handler.h>

template <int dim> class MGDoFCellAccessor;
template <int celldim, int dim> class MGDoFObjectAccessor;

/**
 * Define some types which differ between the dimensions. This class
 * is analogous to the @ref{TriaDimensionInfo} class hierarchy.
 * 
 * @see MGDoFDimensionInfo<1>
 * @see MGDoFDimensionInfo<2>
 */
template <int dim>
class MGDoFDimensionInfo
{
};



/**
 * Define some types for the DoF handling in one dimension.
 *
 * The types have the same meaning as those declared in @ref{TriaDimensionInfo<2>}.
 */
template <>
class MGDoFDimensionInfo<1>
{
  public:
    typedef TriaRawIterator<1,MGDoFCellAccessor<1> >    raw_line_iterator;
    typedef TriaIterator<1,MGDoFCellAccessor<1> >       line_iterator;
    typedef TriaActiveIterator<1,MGDoFCellAccessor<1> > active_line_iterator;

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
 */
template <>
class MGDoFDimensionInfo<2>
{
  public:
    typedef TriaRawIterator<2,MGDoFObjectAccessor<1, 2> >    raw_line_iterator;
    typedef TriaIterator<2,MGDoFObjectAccessor<1, 2> >       line_iterator;
    typedef TriaActiveIterator<2,MGDoFObjectAccessor<1, 2> > active_line_iterator;
    
    typedef TriaRawIterator<2,MGDoFCellAccessor<2> >         raw_quad_iterator;
    typedef TriaIterator<2,MGDoFCellAccessor<2> >            quad_iterator;
    typedef TriaActiveIterator<2,MGDoFCellAccessor<2> >      active_quad_iterator;

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
 * The types have the same meaning as those declared in @ref{TriaDimensionInfo<2>}.
 */
template <>
class MGDoFDimensionInfo<3>
{
  public:
    typedef TriaRawIterator<3,MGDoFObjectAccessor<1, 3> >    raw_line_iterator;
    typedef TriaIterator<3,MGDoFObjectAccessor<1, 3> >       line_iterator;
    typedef TriaActiveIterator<3,MGDoFObjectAccessor<1, 3> > active_line_iterator;

    typedef TriaRawIterator<3,MGDoFObjectAccessor<2, 3> >    raw_quad_iterator;
    typedef TriaIterator<3,MGDoFObjectAccessor<2, 3> >       quad_iterator;
    typedef TriaActiveIterator<3,MGDoFObjectAccessor<2, 3> > active_quad_iterator;

    typedef TriaRawIterator<3,MGDoFCellAccessor<3> >               raw_hex_iterator;
    typedef TriaIterator<3,MGDoFCellAccessor<3> >                  hex_iterator;
    typedef TriaActiveIterator<3,MGDoFCellAccessor<3> >            active_hex_iterator;

    typedef raw_hex_iterator    raw_cell_iterator;
    typedef hex_iterator        cell_iterator;
    typedef active_hex_iterator active_cell_iterator;

    typedef raw_quad_iterator    raw_face_iterator;
    typedef quad_iterator        face_iterator;
    typedef active_quad_iterator active_face_iterator;    
};



/**
 * This class manages degrees of freedom for a multilevel hierarchy of
 * grids. It does mostly the same as does the @p{DoDHandler} class, but
 * it uses a separate enumeration of the degrees of freedom on each
 * level. For example, a vertex has several DoF numbers, one for each
 * level of the triangulation on which it exists.
 *
 * At present, multilevel algorithms are not fully functional, so this
 * documentation is still very brief.
//TODO:[WB] Extend MGDoFHandler doc 
 *
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class MGDoFHandler : public DoFHandler<dim>
{
  public:
    typedef typename MGDoFDimensionInfo<dim>::raw_line_iterator raw_line_iterator;
    typedef typename MGDoFDimensionInfo<dim>::line_iterator line_iterator;
    typedef typename MGDoFDimensionInfo<dim>::active_line_iterator active_line_iterator;

    typedef typename MGDoFDimensionInfo<dim>::raw_quad_iterator raw_quad_iterator;
    typedef typename MGDoFDimensionInfo<dim>::quad_iterator quad_iterator;
    typedef typename MGDoFDimensionInfo<dim>::active_quad_iterator active_quad_iterator;

    typedef typename MGDoFDimensionInfo<dim>::raw_hex_iterator raw_hex_iterator;
    typedef typename MGDoFDimensionInfo<dim>::hex_iterator hex_iterator;
    typedef typename MGDoFDimensionInfo<dim>::active_hex_iterator active_hex_iterator;

    typedef typename MGDoFDimensionInfo<dim>::raw_cell_iterator raw_cell_iterator;
    typedef typename MGDoFDimensionInfo<dim>::cell_iterator cell_iterator;
    typedef typename MGDoFDimensionInfo<dim>::active_cell_iterator active_cell_iterator;

    typedef typename MGDoFDimensionInfo<dim>::raw_face_iterator raw_face_iterator;
    typedef typename MGDoFDimensionInfo<dim>::face_iterator face_iterator;
    typedef typename MGDoFDimensionInfo<dim>::active_face_iterator active_face_iterator;

				     /**
				      * Constructor. Take @p{tria} as the
				      * triangulation to work on.
				      */
    MGDoFHandler (Triangulation<dim> &tria);

				     /**
				      * Destructor
				      */
    virtual ~MGDoFHandler ();
    
				     /**
				      * Go through the triangulation and
				      * distribute the degrees of freedoms
				      * needed for the given finite element
				      * according to the given distribution
				      * method. We first call the @ref{DoFHandler}'s
				      * function and then distribute the
				      * levelwise numbers.
				      *
				      * A copy of the transferred finite
				      * element is stored.
				      *
				      * This function uses the user flags of the
				      * triangulation object, so make sure you
				      * don't need them after calling this
				      * function, or if so store them.
				      */
    virtual void distribute_dofs (const FiniteElement<dim> &, unsigned int offset=0);

				     /**
				      * Actually do the renumbering based on
				      * a list of new dof numbers for all the
				      * dofs. 
				      *
				      * @p{new_numbers} is an array of integers
				      * with size equal to the number of dofs
				      * on the present level. It stores the new
				      * indices after renumbering in the
				      * order of the old indices.
				      */
    void renumber_dofs (const unsigned int               level,
			const std::vector<unsigned int> &new_numbers);

				     /*--------------------------------------*/
    
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
    raw_line_iterator begin_raw_line (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used line
				      *  on level @p{level}.
				      */
    line_iterator     begin_line (const unsigned int level = 0) const;

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
    raw_line_iterator end_line () const;

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
    raw_hex_iterator    begin_raw_hex   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used hex
				      *  on level @p{level}.
				      */
    hex_iterator        begin_hex       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  hex on level @p{level}.
				      */
    active_hex_iterator begin_active_hex(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_hex_iterator    end_hex () const;

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
    raw_hex_iterator    last_raw_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  hex of the level @p{level}, used or not.

				      */
    raw_hex_iterator    last_raw_hex (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used hex.
				      */
    hex_iterator        last_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used hex on level @p{level}.
				      */
    hex_iterator        last_hex (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active hex.
				      */
    active_hex_iterator last_active_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active hex on level @p{level}.
				      */
    active_hex_iterator last_active_hex (const unsigned int level) const;
				     /*@}*/
    
				     /*---------------------------------------*/

    				     /**
				      * Return the number of degrees of freedom
				      * on the specified level.
				      * Included in this number are those
				      * DoFs which are constrained by
				      * hanging nodes.
				      */
    unsigned int n_dofs (const unsigned int level) const;

				     /**
				      * Exception.
				      */
    DeclException1 (ExcInvalidLevel,
		    int,
		    << "The level index " << arg1 << "was not in the valid range.");
    
  private:

				     /** We need for each vertex a list of the
				      * degree of freedom indices on each of the
				      * levels this vertex lives on. Since most
				      * vertices live only on a few levels, it
				      * is not economical to reserve indices for
				      * all the levels there are; rather, we
				      * create an object which holds the indices
				      * on those levels only where the vertex
				      * lives. To construct such an array, it
				      * is necessary to know beforehand which
				      * is the coarsest level the vertex lives
				      * on, how many levels it lives on and
				      * how many dofs there are on each vertex.
				      * If we have this information, we can
				      * allocate exactly the amount of memory
				      * which is needed and need not handle
				      * growing arrays and the like.
				      */
    class MGVertexDoFs {
      public:
	
					 /**
					  * Constructor. This one is empty
					  * because it is difficult to make it
					  * efficient to use vector<>'s and
					  * still construct the object using
					  * the constructor. Use the @p{init}
					  * function to really allocate
					  * memory.
					  */
	MGVertexDoFs ();

					 /**
					  * Allocate memory and
					  * set all indices to @p{-1}.
					  */
	void init (const unsigned int coarsest_level,
		   const unsigned int finest_level,
		   const unsigned int dofs_per_vertex);

					 /**
					  * Destructor
					  */
	~MGVertexDoFs ();

					 /**
					  * Set the index with number
					  * @p{dof_number} of this vertex on
					  * @p{level} to the given index. To
					  * compute the position in the array,
					  * one has to specify how many
					  * dofs per vertex there are. It is
					  * not checked that the level number
					  * is below the number of the finest
					  * level this vertex lives on.
					  *
					  * The function is inline, so should
					  * be reasonably fast.
					  */
	void set_index (const unsigned int level,
			const unsigned int dof_number,
			const unsigned int dofs_per_vertex,
			const unsigned int index);
	
					 /**
					  * Return the index with number
					  * @p{dof_number} of this vertex on
					  * @p{level}. To
					  * compute the position in the array,
					  * one has to specify how many
					  * dofs per vertex there are. It is
					  * not checked that the level number
					  * is below the number of the finest
					  * level this vertex lives on.
					  *
					  * The function is inline, so should
					  * be reasonably fast.
					  */
	unsigned int get_index (const unsigned int level,
				const unsigned int dof_number,
				const unsigned int dofs_per_vertex) const;

					 /**
					  * Return the index of the coarsest
					  * level this vertex lives on.
					  */
	unsigned int get_coarsest_level () const;

					 /**
					  * Return the index of the finest
					  * level this vertex lives on.
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
					 /**
					  * Exception.
					  */
	DeclException0 (ExcInternalError);
	
      private:
					 /**
					  * Store the coarsest level this
					  * vertex lives on. This is used
					  * as an offset when accessing the
					  * dofs of a specific level.
					  */
	unsigned int coarsest_level;

					 /**
					  * Finest level this level lives on.
					  * This is mostly used for error
					  * checking but can also be accessed
					  * by the function @p{get_finest_level}.
					  */
	unsigned int finest_level;

					 /**
					  * Array holding the indices.
					  */
	unsigned int *indices;
    };


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
 *
 * Note that unlike for the usual dofs,
 * here all cells and not only active
 * ones are allowed.
 */
    unsigned int distribute_dofs_on_cell (cell_iterator &cell,
					  unsigned int   next_free_dof);
    
				     /**
				      * Reserve enough space for the MG dof
				      * indices for a given triangulation.
				      */
    void reserve_space ();

    				     /**
				      * Space to store the DoF numbers for the
				      * different levels. Unlike the @p{levels}
				      * object in the @ref{DoFHandler}, these are
				      * not global numbers but rather are
				      * numbers which start from zero on each
				      * level.
				      */
    typename std::vector<DoFLevel<dim>*>    mg_levels;

				     /**
				      * For each vertex there is a list of
				      * indices of the degrees of freedom indices
				      * on the different levels it lives on and
				      * which are these levels.
				      */
    typename std::vector<MGVertexDoFs>      mg_vertex_dofs;
    
				     /**
				      * Vectors storing the number of degrees of
				      * freedom on each level.
				      */
    std::vector<unsigned int>      mg_used_dofs;

#if (__GNUC__==2) && (__GNUC_MINOR__ < 95)  
				     // this seems to be disallowed
				     // by the standard, so gcc2.95
				     // does not accept it
    friend class MGDoFObjectAccessor<1, dim>;
    friend class MGDoFObjectAccessor<2, dim>;
    friend class MGDoFObjectAccessor<3, dim>;
#else
				     // this, however, may grant
				     // access to too many classes,
				     // but ...
    template <int dim1, int dim2> friend class MGDoFObjectAccessor;
#endif
};


/* ----------------------- Inline functions of MGVertexDoFs -------------------*/


template <int dim>
inline
void MGDoFHandler<dim>::MGVertexDoFs::set_index  (const unsigned int level,
						  const unsigned int dof_number,
						  const unsigned int dofs_per_vertex,
						  const unsigned int index) {
  Assert ((level >= coarsest_level) && (level <= finest_level),
	  ExcInvalidLevel(level));
  Assert (dof_number < dofs_per_vertex,
	  ExcIndexRange(dof_number, 0, dofs_per_vertex));
  
  indices[(level-coarsest_level)*dofs_per_vertex + dof_number] = index;
};


template <int dim>
inline
unsigned int
MGDoFHandler<dim>::MGVertexDoFs::get_index  (const unsigned int level,
					     const unsigned int dof_number,
					     const unsigned int dofs_per_vertex) const {
  Assert ((level >= coarsest_level) && (level <= finest_level),
	  ExcInvalidLevel(level));
  Assert (dof_number < dofs_per_vertex,
	  ExcIndexRange (dof_number, 0, dofs_per_vertex));
  
  return indices[(level-coarsest_level)*dofs_per_vertex + dof_number];
};


template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::begin_raw (const unsigned int level) const;
template <>
MGDoFHandler<1>::cell_iterator
MGDoFHandler<1>::begin (const unsigned int level) const;
template <>
MGDoFHandler<1>::active_cell_iterator
MGDoFHandler<1>::begin_active (const unsigned int level) const;
template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::end () const;
template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::last_raw () const;
template <>
MGDoFHandler<1>::raw_cell_iterator
MGDoFHandler<1>::last_raw (const unsigned int level) const;
template <>
MGDoFHandler<1>::cell_iterator
MGDoFHandler<1>::last () const;
template <>
MGDoFHandler<1>::cell_iterator
MGDoFHandler<1>::last (const unsigned int level) const;
template <>
MGDoFHandler<1>::active_cell_iterator
MGDoFHandler<1>::last_active () const;
template <>
MGDoFHandler<1>::active_cell_iterator
MGDoFHandler<1>::last_active (const unsigned int level) const;
template <>
MGDoFDimensionInfo<1>::raw_face_iterator
MGDoFHandler<1>::begin_raw_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<1>::face_iterator
MGDoFHandler<1>::begin_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<1>::active_face_iterator
MGDoFHandler<1>::begin_active_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<1>::raw_face_iterator
MGDoFHandler<1>::end_face () const;
template <>
MGDoFDimensionInfo<1>::raw_face_iterator
MGDoFHandler<1>::last_raw_face () const;
template <>
MGDoFDimensionInfo<1>::raw_face_iterator
MGDoFHandler<1>::last_raw_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<1>::face_iterator
MGDoFHandler<1>::last_face () const;
template <>
MGDoFDimensionInfo<1>::face_iterator
MGDoFHandler<1>::last_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<1>::active_face_iterator
MGDoFHandler<1>::last_active_face () const;
template <>
MGDoFDimensionInfo<1>::active_face_iterator
MGDoFHandler<1>::last_active_face (const unsigned int) const;
template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::begin_raw_quad (const unsigned int) const;
template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::begin_quad (const unsigned int) const;
template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::begin_active_quad (const unsigned int) const;
template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::end_quad () const;
template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::last_raw_quad (const unsigned int) const;
template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::last_quad (const unsigned int) const;
template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::last_active_quad (const unsigned int) const;
template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::last_raw_quad () const;
template <>
MGDoFHandler<1>::raw_quad_iterator
MGDoFHandler<1>::last_raw_quad () const;
template <>
MGDoFHandler<1>::quad_iterator
MGDoFHandler<1>::last_quad () const;
template <>
MGDoFHandler<1>::active_quad_iterator
MGDoFHandler<1>::last_active_quad () const;
template <>
MGDoFHandler<1>::raw_hex_iterator
MGDoFHandler<1>::begin_raw_hex (const unsigned int) const;
template <>
MGDoFHandler<1>::hex_iterator
MGDoFHandler<1>::begin_hex (const unsigned int) const;
template <>
MGDoFHandler<1>::active_hex_iterator
MGDoFHandler<1>::begin_active_hex (const unsigned int) const;
template <>
MGDoFHandler<1>::raw_hex_iterator
MGDoFHandler<1>::end_hex () const;
template <>
MGDoFHandler<1>::raw_hex_iterator
MGDoFHandler<1>::last_raw_hex (const unsigned int) const;
template <>
MGDoFHandler<1>::hex_iterator
MGDoFHandler<1>::last_hex (const unsigned int) const;
template <>
MGDoFHandler<1>::active_hex_iterator
MGDoFHandler<1>::last_active_hex (const unsigned int) const;
template <>
MGDoFHandler<1>::raw_hex_iterator
MGDoFHandler<1>::last_raw_hex () const;
template <>
MGDoFHandler<1>::hex_iterator
MGDoFHandler<1>::last_hex () const;
template <>
MGDoFHandler<1>::active_hex_iterator
MGDoFHandler<1>::last_active_hex () const;
//======================================================================//
template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::begin_raw (const unsigned int level) const;
template <>
MGDoFHandler<2>::cell_iterator
MGDoFHandler<2>::begin (const unsigned int level) const;
template <>
MGDoFHandler<2>::active_cell_iterator
MGDoFHandler<2>::begin_active (const unsigned int level) const;
template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::end () const;
template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::last_raw () const;
template <>
MGDoFHandler<2>::raw_cell_iterator
MGDoFHandler<2>::last_raw (const unsigned int level) const;
template <>
MGDoFHandler<2>::cell_iterator
MGDoFHandler<2>::last () const;
template <>
MGDoFHandler<2>::cell_iterator
MGDoFHandler<2>::last (const unsigned int level) const;
template <>
MGDoFHandler<2>::active_cell_iterator
MGDoFHandler<2>::last_active () const;
template <>
MGDoFHandler<2>::active_cell_iterator
MGDoFHandler<2>::last_active (const unsigned int level) const;
template <>
MGDoFDimensionInfo<2>::raw_face_iterator
MGDoFHandler<2>::begin_raw_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<2>::face_iterator
MGDoFHandler<2>::begin_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<2>::active_face_iterator
MGDoFHandler<2>::begin_active_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<2>::raw_face_iterator
MGDoFHandler<2>::end_face () const;
template <>
MGDoFDimensionInfo<2>::raw_face_iterator
MGDoFHandler<2>::last_raw_face () const;
template <>
MGDoFDimensionInfo<2>::raw_face_iterator
MGDoFHandler<2>::last_raw_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<2>::face_iterator
MGDoFHandler<2>::last_face () const;
template <>
MGDoFDimensionInfo<2>::face_iterator
MGDoFHandler<2>::last_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<2>::active_face_iterator
MGDoFHandler<2>::last_active_face () const;
template <>
MGDoFDimensionInfo<2>::active_face_iterator
MGDoFHandler<2>::last_active_face (const unsigned int) const;
template <>
MGDoFHandler<2>::raw_hex_iterator
MGDoFHandler<2>::begin_raw_hex (const unsigned int) const;
template <>
MGDoFHandler<2>::hex_iterator
MGDoFHandler<2>::begin_hex (const unsigned int) const;
template <>
MGDoFHandler<2>::active_hex_iterator
MGDoFHandler<2>::begin_active_hex (const unsigned int) const;
template <>
MGDoFHandler<2>::raw_hex_iterator
MGDoFHandler<2>::end_hex () const;
template <>
MGDoFHandler<2>::raw_hex_iterator
MGDoFHandler<2>::last_raw_hex (const unsigned int) const;
template <>
MGDoFHandler<2>::hex_iterator
MGDoFHandler<2>::last_hex (const unsigned int) const;
template <>
MGDoFHandler<2>::active_hex_iterator
MGDoFHandler<2>::last_active_hex (const unsigned int) const;
template <>
MGDoFHandler<2>::raw_hex_iterator
MGDoFHandler<2>::last_raw_hex () const;
template <>
MGDoFHandler<2>::hex_iterator
MGDoFHandler<2>::last_hex () const;
template <>
MGDoFHandler<2>::active_hex_iterator
MGDoFHandler<2>::last_active_hex () const;
//======================================================================//
template <>
MGDoFHandler<3>::raw_cell_iterator
MGDoFHandler<3>::begin_raw (const unsigned int level) const;
template <>
MGDoFHandler<3>::cell_iterator
MGDoFHandler<3>::begin (const unsigned int level) const;
template <>
MGDoFHandler<3>::active_cell_iterator
MGDoFHandler<3>::begin_active (const unsigned int level) const;
template <>
MGDoFHandler<3>::raw_cell_iterator
MGDoFHandler<3>::end () const;
template <>
MGDoFHandler<3>::raw_cell_iterator
MGDoFHandler<3>::last_raw () const;
template <>
MGDoFHandler<3>::raw_cell_iterator
MGDoFHandler<3>::last_raw (const unsigned int level) const;
template <>
MGDoFHandler<3>::cell_iterator
MGDoFHandler<3>::last () const;
template <>
MGDoFHandler<3>::cell_iterator
MGDoFHandler<3>::last (const unsigned int level) const;
template <>
MGDoFHandler<3>::active_cell_iterator
MGDoFHandler<3>::last_active () const;
template <>
MGDoFHandler<3>::active_cell_iterator
MGDoFHandler<3>::last_active (const unsigned int level) const;
template <>
MGDoFDimensionInfo<3>::raw_face_iterator
MGDoFHandler<3>::begin_raw_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<3>::face_iterator
MGDoFHandler<3>::begin_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<3>::active_face_iterator
MGDoFHandler<3>::begin_active_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<3>::raw_face_iterator
MGDoFHandler<3>::end_face () const;
template <>
MGDoFDimensionInfo<3>::raw_face_iterator
MGDoFHandler<3>::last_raw_face () const;
template <>
MGDoFDimensionInfo<3>::raw_face_iterator
MGDoFHandler<3>::last_raw_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<3>::face_iterator
MGDoFHandler<3>::last_face () const;
template <>
MGDoFDimensionInfo<3>::face_iterator
MGDoFHandler<3>::last_face (const unsigned int) const;
template <>
MGDoFDimensionInfo<3>::active_face_iterator
MGDoFHandler<3>::last_active_face () const;
template <>
MGDoFDimensionInfo<3>::active_face_iterator
MGDoFHandler<3>::last_active_face (const unsigned int) const;


template <>
void MGDoFHandler<1>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers);
template <>
void MGDoFHandler<2>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers);
template <>
void MGDoFHandler<3>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers);


/*----------------------------   mg_dof.h     ---------------------------*/

#endif
/*----------------------------   mg_dof.h     ---------------------------*/
