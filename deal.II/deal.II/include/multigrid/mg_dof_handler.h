/*----------------------------   mg_dof.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __mg_dof_H
#define __mg_dof_H
/*----------------------------   mg_dof.h     ---------------------------*/


#include <grid/dof.h>



/**
 * Define some types which differ between the dimensions. This class
 * is analogous to the \Ref{TriaDimensionInfo} class hierarchy.
 * 
 * @see MGDoFDimensionInfo<1>
 * @see MGDoFDimensionInfo<2>
 */
template <int dim>
class MGDoFDimensionInfo;





/**
 * Define some types for the DoF handling in one dimension.
 *
 * The types have the same meaning as those declared in \Ref{TriaDimensionInfo<2>}.
 */
class MGDoFDimensionInfo<1> {
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
 * The types have the same meaning as those declared in \Ref{TriaDimensionInfo<2>}.
 */
class MGDoFDimensionInfo<2> {
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
 * The types have the same meaning as those declared in \Ref{TriaDimensionInfo<2>}.
 */
class MGDoFDimensionInfo<3> {
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





template <int dim>
class MGDoFHandler
  :
  public DoFHandler<dim>
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
				      * Constructor. Take #tria# as the
				      * triangulation to work on.
				      */
    MGDoFHandler (Triangulation<dim> *tria);

				     /**
				      * Destructor
				      */
    virtual ~MGDoFHandler ();
    
				     /**
				      * Go through the triangulation and
				      * distribute the degrees of freedoms
				      * needed for the given finite element
				      * according to the given distribution
				      * method. We first call the #DoFHandler#'s
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
				      * #new_numbers# is an array of integers
				      * with size equal to the number of dofs
				      * on the present level. It stores the new
				      * indices after renumbering in the
				      * order of the old indices.
				      */
    void renumber_dofs (const unsigned int level,
			const vector<int> &new_numbers);

				     /**
				      * Write the sparsity structure of the
				      * matrix belonging to the specified
				      * #level# including constrained
				      * degrees of freedom into the
				      * matrix structure. The sparsity pattern
				      * does not include entries introduced by
				      * the elimination of constrained nodes.
				      * The sparsity
				      * pattern is not compressed, since if
				      * you want to call
				      * #ConstraintMatrix::condense(1)#
				      * afterwards, new entries have to be
				      * added. However, if you don't want to call
				      * #ConstraintMatrix::condense(1)#, you
				      * have to compress the matrix yourself,
				      * using #SparseMatrixStruct::compress()#.
				      */
    void make_sparsity_pattern (const unsigned int  level,
				SparseMatrixStruct &sparsity) const; 

				     /*--------------------------------------*/
    
				     /**
				      *  @name Cell iterator functions
				      */
				     /*@{*/
				     /**
				      *  Iterator to the first cell, used
				      *  or not, on level #level#. If a level
				      *  has no cells, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls #begin_raw_line#
				      *  in 1D and #begin_raw_quad# in 2D.
				      */
    raw_cell_iterator    begin_raw   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used cell
				      *  on level #level#.
				      *
				      *  This function calls #begin_line#
				      *  in 1D and #begin_quad# in 2D.
				      */
    cell_iterator        begin       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  cell on level #level#.
				      *
				      *  This function calls #begin_active_line#
				      *  in 1D and #begin_active_quad# in 2D.
				      */
    active_cell_iterator begin_active(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls #end_line#
				      *  in 1D and #end_quad# in 2D.
				      */
    raw_cell_iterator    end () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    cell_iterator        end (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_cell_iterator    end_raw (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_cell_iterator end_active (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last cell, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 1D and #last_raw_quad# in 2D.
				      */
    raw_cell_iterator    last_raw () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  cell of the level #level#, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 1D and #last_raw_quad# in 2D.
				      */
    raw_cell_iterator    last_raw (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell.
				      *
				      *  This function calls #last_line#
				      *  in 1D and #last_quad# in 2D.
				      */
    cell_iterator        last () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell on level #level#.
				      *
				      *  This function calls #last_line#
				      *  in 1D and #last_quad# in 2D.
				      */
    cell_iterator        last (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active cell.
				      *
				      *  This function calls #last_active_line#
				      *  in 1D and #last_active_quad# in 2D.
				      */
    active_cell_iterator last_active () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active cell on level #level#.
				      *
				      *  This function calls #last_active_line#
				      *  in 1D and #last_active_quad# in 2D.
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
				      *  or not, on level #level#. If a level
				      *  has no faces, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls #begin_raw_line#
				      *  in 2D and #begin_raw_quad# in 3D.
				      */
    raw_face_iterator    begin_raw_face   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used face
				      *  on level #level#.
				      *
				      *  This function calls #begin_line#
				      *  in 2D and #begin_quad# in 3D.
				      */
    face_iterator        begin_face       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  face on level #level#.
				      *
				      *  This function calls #begin_active_line#
				      *  in 2D and #begin_active_quad# in 3D.
				      */
    active_face_iterator begin_active_face(const unsigned int level = 0) const;

				     /**
				      *  Iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls #end_line#
				      *  in 2D and #end_quad# in 3D.
				      */
    raw_face_iterator    end_face () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    face_iterator        end_face (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_face_iterator    end_raw_face (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_face_iterator end_active_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last face, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 2D and #last_raw_quad# in 3D.
				      */
    raw_face_iterator    last_raw_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  face of the level #level#, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 2D and #last_raw_quad# in 3D.
				      */
    raw_face_iterator    last_raw_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face.
				      *
				      *  This function calls #last_line#
				      *  in 2D and #last_quad# in 3D.
				      */
    face_iterator        last_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face on level #level#.
				      *
				      *  This function calls #last_line#
				      *  in 2D and #last_quad# in 3D.
				      */
    face_iterator        last_face (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active face.
				      *
				      *  This function calls #last_active_line#
				      *  in 2D and #last_active_quad# in 3D.
				      */
    active_face_iterator last_active_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active face on level #level#.
				      *
				      *  This function calls #last_active_line#
				      *  in 2D and #last_active_quad# in 3D.
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
				      *  or not, on level #level#. If a level
				      *  has no lines, a past-the-end iterator
				      *  is returned.
				      */
    raw_line_iterator begin_raw_line (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used line
				      *  on level #level#.
				      */
    line_iterator     begin_line (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  line on level #level#.
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
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    line_iterator        end_line (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_line_iterator    end_raw_line (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_line_iterator end_active_line (const unsigned int level) const;


				     /**
				      *  Return an iterator pointing to the
				      *  last line, used or not.
				      */
    raw_line_iterator    last_raw_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  line of the level #level#, used or not.

				      */
    raw_line_iterator    last_raw_line (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line.
				      */
    line_iterator        last_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line on level #level#.
				      */
    line_iterator        last_line (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active line.
				      */
    active_line_iterator last_active_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active line on level #level#.
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
				      *  or not, on level #level#. If a level
				      *  has no quads, a past-the-end iterator
				      *  is returned.
				      */
    raw_quad_iterator    begin_raw_quad   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used quad
				      *  on level #level#.
				      */
    quad_iterator        begin_quad       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  quad on level #level#.
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
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    quad_iterator        end_quad (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_quad_iterator    end_raw_quad (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_quad_iterator end_active_quad (const unsigned int level) const;


				     /**
				      *  Return an iterator pointing to the
				      *  last quad, used or not.
				      */
    raw_quad_iterator    last_raw_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  quad of the level #level#, used or not.

				      */
    raw_quad_iterator    last_raw_quad (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad.
				      */
    quad_iterator        last_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad on level #level#.
				      */
    quad_iterator        last_quad (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active quad.
				      */
    active_quad_iterator last_active_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active quad on level #level#.
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
				      *  or not, on level #level#. If a level
				      *  has no hexs, a past-the-end iterator
				      *  is returned.
				      */
    raw_hex_iterator    begin_raw_hex   (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first used hex
				      *  on level #level#.
				      */
    hex_iterator        begin_hex       (const unsigned int level = 0) const;

				     /**
				      *  Iterator to the first active
				      *  hex on level #level#.
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
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    hex_iterator        end_hex (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_hex_iterator    end_raw_hex (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_hex_iterator end_active_hex (const unsigned int level) const;


				     /**
				      *  Return an iterator pointing to the
				      *  last hex, used or not.
				      */
    raw_hex_iterator    last_raw_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  hex of the level #level#, used or not.

				      */
    raw_hex_iterator    last_raw_hex (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used hex.
				      */
    hex_iterator        last_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used hex on level #level#.
				      */
    hex_iterator        last_hex (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active hex.
				      */
    active_hex_iterator last_active_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active hex on level #level#.
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
					  * the constructor. Use the #init#
					  * function to really allocate
					  * memory.
					  */
	MGVertexDoFs ();

					 /**
					  * Allocate memory and
					  * set all indices to #-1#.
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
					  * #dof_number# of this vertex on
					  * #level# to the given index. To
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
					  * #dof_number# of this vertex on
					  * #level#. To
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
	int get_index (const unsigned int level,
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
					  * by the function #get_finest_level#.
					  */
	unsigned int finest_level;

					 /**
					  * Array holding the indices.
					  */
	int         *indices;
    };


				     /**
				      * Distribute dofs on the given cell,
				      * with new dofs starting with index
				      * #next_free_dof#. Return the next
				      * unused index number. The finite
				      * element used is the one given to
				      * #distribute_dofs#, which is copied
				      * to #selected_fe#.
				      *
				      * This function is excluded from the
				      * #distribute_dofs# function since
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
				      * different levels. Unlike the #levels#
				      * object in the #DoFHandler#, these are
				      * not global numbers but rather are
				      * numbers which start from zero on each
				      * level.
				      */
    vector<DoFLevel<dim>*>    mg_levels;

				     /**
				      * For each vertex there is a list of
				      * indices of the degrees of freedom indices
				      * on the different levels it lives on and
				      * which are these levels.
				      */
    vector<MGVertexDoFs>      mg_vertex_dofs;
    
				     /**
				      * Vectors storing the number of degrees of
				      * freedom on each level.
				      */
    vector<unsigned int>      mg_used_dofs;

    friend class MGDoFObjectAccessor<1, dim>;
    friend class MGDoFObjectAccessor<2, dim>;
    friend class MGDoFObjectAccessor<3, dim>;
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
int MGDoFHandler<dim>::MGVertexDoFs::get_index  (const unsigned int level,
						 const unsigned int dof_number,
						 const unsigned int dofs_per_vertex) const {
  Assert ((level >= coarsest_level) && (level <= finest_level),
	  ExcInvalidLevel(level));
  Assert (dof_number < dofs_per_vertex,
	  ExcIndexRange (dof_number, 0, dofs_per_vertex));
  
  return indices[(level-coarsest_level)*dofs_per_vertex + dof_number];
};



/*----------------------------   mg_dof.h     ---------------------------*/
/* end of #ifndef __mg_dof_H */
#endif
/*----------------------------   mg_dof.h     ---------------------------*/
