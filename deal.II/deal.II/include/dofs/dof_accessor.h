/*----------------------------   dof_iterator.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __dof_iterator_H
#define __dof_iterator_H
/*----------------------------   dof_iterator.h     ---------------------------*/


#include <grid/tria_accessor.h>
#include <lac/forward-declarations.h>
#include <vector>


// note: in non-debug mode, i.e. with optimizations, the file
// dof_accessor.templates.h is included at the end of this file.
// this includes a lot of templates and thus makes compilation
// slower, but at the same time allows for more aggressive
// inlining and thus faster code.



/**
 * Define the basis for accessors to the degrees of freedom.
 *
 * Note that it is allowed to construct an object of which the
 * #dof_handler# pointer is a Null pointer. Such an object would
 * result in a strange kind of behaviour, though every reasonable
 * operating system should disallow access through that pointer.
 * The reason we do not check for the null pointer in the
 * constructor which gets passed the #DoFHandler# pointer is that
 * if we did we could not make dof iterators member of other classes
 * (like in the #FEValues# class) if we did not know about the
 * #DoFHandler# object to be used upon construction of that object.
 * Through the way this class is implemented here, we allow the
 * creation of a kind of virgin object which only gets useful if
 * assigned to from another object before first usage.
 *
 * Opposite to construction, it is not possible to copy an object
 * which has an invalid dof handler pointer. This is to guarantee
 * that every iterator which is once assigned to is a valid
 * object. However, this assertion only holds in debug mode, when
 * the #Assert# macro is switched on.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class DoFAccessor {
  public:
				     /**
				      * Constructor
				      */
    DoFAccessor () : dof_handler(0) {
      Assert (false, ExcInvalidObject());
    };

				     /**
				      * This should be the default constructor.
				      * We cast away the #const#ness of the
				      * pointer which clearly is EVIL but
				      * we can't help without making all
				      * functions which could somehow use
				      * iterators (directly or indirectly) make
				      * non-const, even if they preserve
				      * constness.
				      */
    DoFAccessor (const DoFHandler<dim> *dof_handler) :
		    dof_handler(const_cast<DoFHandler<dim>*>(dof_handler)) {};

				     /**
				      * Reset the DoF handler pointer.
				      */
    void set_dof_handler (DoFHandler<dim> *dh) {
      Assert (dh != 0, ExcInvalidObject());
      dof_handler = dh;
    };

				     /**
				      * Copy operator.
				      */
    DoFAccessor<dim> & operator = (const DoFAccessor<dim> &da) {
      set_dof_handler (da.dof_handler);
      return *this;
    };
    
				     /**
				      * Exception for child classes
				      */
    DeclException0 (ExcInvalidObject);
				     /**
				      * Exception
				      */
    DeclException0 (ExcVectorNotEmpty);
				     /**
				      * Exception
				      */
    DeclException0 (ExcVectorDoesNotMatch);
				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixDoesNotMatch);

  protected:
				     /**
				      * Store the address of the #DoFHandler# object
				      * to be accessed.
				      */
    DoFHandler<dim> *dof_handler;  
};



/**
 * Common template for line, quad, hex.
 *
 * Internal: inheritance is necessary for the general template due to
 * a compiler error.
 * @author Guido Kanschat, 1999
 */
template<int celldim, int dim, typename BaseClass>
class DoFObjectAccessor : public DoFAccessor<dim>,
			  public BaseClass
{};


/**
 * Closure class.
 */
template<int dim, typename BaseClass>
class DoFObjectAccessor<0, dim, BaseClass> : public DoFAccessor<dim>,
					     public BaseClass
{
  public:
    typedef void* AccessorData;
    DoFObjectAccessor (Triangulation<dim> *,
		     const int,
		     const int,
		     const AccessorData *)
      {}
};


/**
 * Grant access to the degrees of freedom located on lines.
 * This class follows mainly the route laid out by the accessor library
 * declared in the triangulation library (\Ref{TriaAccessor}). It enables
 * the user to access the degrees of freedom on the lines (there are similar
 * versions for the DoFs on quads, etc), where the dimension of the underlying
 * triangulation does not really matter (i.e. this accessor works with the
 * lines in 1D-, 2D-, etc dimensions).
 *
 *
 * \subsection{Usage}
 *
 * The \Ref{DoFDimensionInfo} classes inherited by the \Ref{DoFHandler} classes
 * declare typedefs to iterators using the accessors declared in this class
 * hierarchy tree. Usage is best to happens through these typedefs, since they
 * are more secure to changes in the class naming and template interface as well
 * as they provide easier typing (much less complicated names!).
 * 
 * 
 * \subsection{Notes about the class hierarchy structure}
 *
 * The class hierarchy seems to be a bit confused here. The reason for this is
 * that we would really like to derive a #DoFLineAccessor# from a #LineAccessor#.
 * Unfortunately, we would run into problems, if we wanted a #DoFLineAccessor#
 * in one spatial dimension, in which case a line is also a cell. The traditional
 * solution would be to declare a #DoFCellAccessor<1># which is derived from
 * #DoFObjectAccessor<1, 1># and #CellAccessor<1># (the #DoFObjectAccessor<1, dim># cannot
 * itself be derived from #CellAccessor<dim># since a line is not a cell
 * unless in one space dimension), but since a #DoFLineAccessor# and a
 * #CellAccessor# are both derived from #TriaAccessor#, we would have to make
 * the last derivation virtual.
 *
 * Since we want to avoid virtual inheritance since this involves another
 * indirection in every member variable access, we chose another way: we
 * pass a second template parameter to a #DoFLineAccessor# which tells it
 * which class to be derived from: if we are in one spatial dimension, the
 * base class is to be #CellAccessor<1>#, in two or more dimensions it
 * is a #TriaObjectAccessor<1, dim>#, i.e. am accessor to lines without the missing
 * functionality needed for cells (neighbors, etc.).
 *
 * This way we can declare a #DoFCellAccessor# in one dimension by deriving
 * from #DoFObjectAccessor<1, 1,CellAccessor<1> >#, thus getting the cell
 * functionality through the #DoFLineAccessor# instead of through a virtual
 * multiple inheritance of #DoFLineAccessor# and #CellAccessor<1>#.
 *
 * The same concept is used with #DoFQuadAccessor# classes etc.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename BaseClass>
class DoFObjectAccessor<1, dim, BaseClass> :  public DoFAccessor<dim>,
					      public BaseClass
{
  public:
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef DoFHandler<dim> AccessorData;
    
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFObjectAccessor ();
    
    				     /**
				      * Constructor. The #local_data#
				      * argument is assumed to be a pointer
				      * to a #DoFHandler<dim># object.
				      */
    DoFObjectAccessor (Triangulation<dim> *tria,
		     const int           level,
		     const int           index,
		     const AccessorData *local_data) :
		    DoFAccessor<dim> (local_data),
		    BaseClass(tria,level,index) {};
    
				     /**
				      * Return the index of the #i#th degree
				      * of freedom of this line.
				      */
    int dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the #i#th degree
				      * of freedom of this line to #index#.
				      */
    void set_dof_index (const unsigned int i, const int index) const;

				     /**
				      * Return the index of the #i#th degree
				      * on the #vertex#th vertex.
				      */
    int vertex_dof_index (const unsigned int vertex,
			  const unsigned int i) const;

				     /**
				      * Set the index of the #i#th degree
				      * on the #vertex#th vertex to #index#.
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const int          index) const;

    				     /**
				      * Return the indices of the dofs of this
				      * line in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, 
				      * dofs on line.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand.
				      */
    void get_dof_indices (vector<int> &dof_indices) const;

				     /**
				      * Return the #i#th child as a DoF line
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a line accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,DoFObjectAccessor<1, dim,BaseClass> > child (const unsigned int) const;

				     /**
				      * Distribute a local (cell based) vector
				      * to a global one by mapping the local
				      * numbering of the degrees of freedom
				      * to the global one and entering the
				      * local values into the global vector.
				      *
				      * The elements are {\it added} up to
				      * the elements in the global vector,
				      * rather than just set, since this is
				      * usually what one wants.
				      */
    void distribute_local_to_global (const Vector<double> &local_source,
				     Vector<double>       &global_destination) const;

				     /**
				      * This function does much the same as the
				      * #distribute_local_to_global(dVector,dVector)#
				      * function, but operates on matrices
				      * instead of vectors. The sparse matrix
				      * is supposed to have non-zero entry
				      * slots where required.
				      */
    void distribute_local_to_global (const FullMatrix<double> &local_source,
				     SparseMatrix<double>     &global_destination) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFObjectAccessor<1, dim,BaseClass> &a);
};



/**
 * Grant access to the degrees of freedom located on quads.
 *
 * @see DoFLineAccessor
 */
template <int dim, typename BaseClass>
class DoFObjectAccessor<2, dim, BaseClass> :  public DoFAccessor<dim>,
					      public BaseClass
{
  public:
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef DoFHandler<dim> AccessorData;
    
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFObjectAccessor ();

    				     /**
				      * Constructor. The #local_data#
				      * argument is assumed to be a pointer
				      * to a #DoFHandler<dim># object.
				      */
    DoFObjectAccessor (Triangulation<dim> *tria,
		     const int           level,
		     const int           index,
		     const AccessorData *local_data) :
		    DoFAccessor<dim> (local_data),
		    BaseClass        (tria,level,index) {};
    
				     /**
				      * Return the index of the #i#th degree
				      * of freedom of this quad.
				      */
    int dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the #i#th degree
				      * of freedom of this quad to #index#.
				      */
    void set_dof_index (const unsigned int i, const int index) const;

				     /**
				      * Return the index of the #i#th degree
				      * on the #vertex#th vertex.
				      */
    int vertex_dof_index (const unsigned int vertex,
			  const unsigned int i) const;

				     /**
				      * Set the index of the #i#th degree
				      * on the #vertex#th vertex to #index#.
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const int          index) const;

    				     /**
				      * Return the indices of the dofs of this
				      * quad in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand.
				      */
    void get_dof_indices (vector<int> &dof_indices) const;

    				     /**
				      *  Return a pointer to the #i#th line
				      *  bounding this #Quad#.
				      */
    TriaIterator<dim,DoFObjectAccessor<1, dim,TriaObjectAccessor<1, dim> > >
    line (const unsigned int i) const;
    
				     /**
				      * Return the #i#th child as a DoF quad
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a quad accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,DoFObjectAccessor<2, dim, BaseClass> > child (const unsigned int) const;

				     /**
				      * Distribute a local (cell based) vector
				      * to a global one by mapping the local
				      * numbering of the degrees of freedom
				      * to the global one and entering the
				      * local values into the global vector.
				      *
				      * The elements are {\it added} up to
				      * the elements in the global vector,
				      * rather than just set, since this is
				      * usually what one wants.
				      */
    void distribute_local_to_global (const Vector<double> &local_source,
				     Vector<double>       &global_destination) const;

				     /**
				      * This function does much the same as the
				      * #distribute_local_to_global(dVector,dVector)#
				      * function, but operates on matrices
				      * instead of vectors. The sparse matrix
				      * is supposed to have non-zero entry
				      * slots where required.
				      */
    void distribute_local_to_global (const FullMatrix<double> &local_source,
				     SparseMatrix<double>     &global_destination) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFObjectAccessor<2, dim,BaseClass> &a);
};





/**
 * Grant access to the degrees of freedom located on quads.
 *
 * @see DoFLineAccessor
 */
template <int dim, typename BaseClass>
class DoFObjectAccessor<3, dim, BaseClass> :  public DoFAccessor<dim>,
					      public BaseClass
{
  public:
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef DoFHandler<dim> AccessorData;
    
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFObjectAccessor ();

    				     /**
				      * Constructor. The #local_data#
				      * argument is assumed to be a pointer
				      * to a #DoFHandler<dim># object.
				      */
    DoFObjectAccessor (Triangulation<dim> *tria,
		    const int           level,
		    const int           index,
		    const AccessorData *local_data) :
		    DoFAccessor<dim> (local_data),
		    BaseClass        (tria,level,index) {};
    
				     /**
				      * Return the index of the #i#th degree
				      * of freedom of this quad.
				      */
    int dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the #i#th degree
				      * of freedom of this quad to #index#.
				      */
    void set_dof_index (const unsigned int i, const int index) const;

				     /**
				      * Return the index of the #i#th degree
				      * on the #vertex#th vertex.
				      */
    int vertex_dof_index (const unsigned int vertex,
			  const unsigned int i) const;

				     /**
				      * Set the index of the #i#th degree
				      * on the #vertex#th vertex to #index#.
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const int          index) const;

    				     /**
				      * Return the indices of the dofs of this
				      * quad in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand.
				      */
    void get_dof_indices (vector<int> &dof_indices) const;

    				     /**
				      *  Return a pointer to the #i#th line
				      *  bounding this #Hex#.
				      */
    TriaIterator<dim,DoFObjectAccessor<1, dim,TriaObjectAccessor<1, dim> > >
    line (const unsigned int i) const;

    				     /**
				      *  Return a pointer to the #i#th quad
				      *  bounding this #Hex#.
				      */
    TriaIterator<dim,DoFObjectAccessor<2, dim,TriaObjectAccessor<2, dim> > >
    quad (const unsigned int i) const;
    
				     /**
				      * Return the #i#th child as a DoF hex
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a hex accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,DoFObjectAccessor<3, dim,BaseClass> > child (const unsigned int) const;

				     /**
				      * Distribute a local (cell based) vector
				      * to a global one by mapping the local
				      * numbering of the degrees of freedom
				      * to the global one and entering the
				      * local values into the global vector.
				      *
				      * The elements are {\it added} up to
				      * the elements in the global vector,
				      * rather than just set, since this is
				      * usually what one wants.
				      */
    void distribute_local_to_global (const Vector<double> &local_source,
				     Vector<double>       &global_destination) const;

				     /**
				      * This function does much the same as the
				      * #distribute_local_to_global(dVector,dVector)#
				      * function, but operates on matrices
				      * instead of vectors. The sparse matrix
				      * is supposed to have non-zero entry
				      * slots where required.
				      */
    void distribute_local_to_global (const FullMatrix<double> &local_source,
				     SparseMatrix<double>     &global_destination) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFObjectAccessor<3, dim,BaseClass> &a);
};



/**
 * Grant access to the degrees of freedom on a cell. In fact, since all
 * access to the degrees of freedom has been enabled by the classes
 * #DoFObjectAccessor<1, 1># and #DoFObjectAccessor<2, 2># for the space dimension
 * one and two, respectively, this class only collects the pieces
 * together by deriving from the appropriate #DoF*Accessor# and the
 * right #CellAccessor<dim># and finally adding two functions which give
 * access to the neighbors and children as #DoFCellAccessor# objects
 * rather than #CellAccessor# objects (the latter function was inherited
 * from the #CellAccessor<dim># class).
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class DoFCellAccessor :  public DoFObjectAccessor<dim, dim, CellAccessor<dim> > {
  public:
    				     /**
				      * Type of faces.
				      */
    typedef
    TriaIterator<dim, DoFObjectAccessor<dim-1, dim,TriaObjectAccessor<dim-1, dim> > >
    face_iterator;

				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef typename DoFObjectAccessor<dim, dim, CellAccessor<dim> >::AccessorData AccessorData;
    
    				     /**
				      * Constructor
				      */
    DoFCellAccessor (Triangulation<dim> *tria,
		     const int           level,
		     const int           index,
		     const AccessorData *local_data) :
		     DoFObjectAccessor<dim, dim, CellAccessor<dim> > (tria,level,index,local_data) {};

				     /**
				      * Return the #i#th neighbor as a DoF cell
				      * iterator. This function is needed since
				      * the neighbor function of the base
				      * class returns a cell accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,DoFCellAccessor<dim> > neighbor (const unsigned int) const;

    				     /**
				      * Return the #i#th child as a DoF cell
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a cell accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,DoFCellAccessor<dim> > child (const unsigned int) const;

    				     /**
				      * Return an iterator to the #i#th face
				      * of this cell.
				      *
				      * This function is not implemented in 1D,
				      * and maps to DoFObjectAccessor<2, dim>::line in 2D.
				      */
    face_iterator
    face (const unsigned int i) const;

    				     /**
				      * Return the values of the given vector
				      * restricted to the dofs of this
				      * cell in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand. This
				      * function is only callable for active
				      * cells.
				      */
    template <typename number>
    void get_dof_values (const Vector<number> &values,
			 Vector<number>       &local_values) const;

				     /**
				      * Return the interpolation of the given
				      * finite element function to the present
				      * cell. In the simplest case, the cell
				      * is a terminal one, i.e. has no children;
				      * then, the returned value is the vector
				      * of nodal values on that cell. You could
				      * then as well get the desired values
				      * through the #get_dof_values# function
				      * In the other case, when the cell has
				      * children, we use the restriction
				      * matrices provided by the finite element
				      * class to compute the interpolation
				      * from the children to the present cell.
				      *
				      * It is assumed that both vectors already
				      * have the right size beforehand. This
				      * function assumes the existence of an
				      * interpolation from child cells to the
				      * mother cell, denoted by the restriction
				      * matrices of the finite element class.
				      * Futhermore, this interpolation should
				      * be possible for each child alone, i.e.
				      * it should be possible to compute the
				      * restriction by writing values obtained
				      * from each child directly into the output
				      * vector, without, for example, computing
				      * an average over all children.
				      * These properties, however, do not exist
				      * for all elements; an example is the
				      * DG(0) element, which represents
				      * piecewise constant elements: for these
				      * the restriction to mother cell could
				      * be the average of the values of the
				      * children, maybe weighted by the measure
				      * of each child. It is not yet decided
				      * what the this function does in these
				      * cases.
				      */
    template <typename number>
    void get_interpolated_dof_values (const Vector<number> &values,
				      Vector<number>       &interpolated_values) const;

				     /**
				      * This function is the counterpart to
				      * #get_dof_values#: it takes a vector
				      * of values for the degrees of freedom
				      * of the cell pointed to by this iterator
				      * and writes these values into the global
				      * data vector #values#. This function
				      * is only callable for active cells.
				      *
				      * Note that for continuous finite
				      * elements, calling this function affects
				      * the dof values on neighboring cells as
				      * well. It may also violate continuity
				      * requirements for hanging nodes, if
				      * neighboring cells are less refined than
				      * the present one. These requirements
				      * are not taken care of and must be
				      * enforced by the user afterwards.
				      *
				      * It is assumed that both vectors already
				      * have the right size beforehand.
				      */
    template <typename number>
    void set_dof_values (const Vector<number> &local_values,
			 Vector<number>       &values) const;

				     /**
				      * This, again, is the counterpart to
				      * #get_interpolated_dof_values#: you
				      * specify the dof values on a cell and
				      * these are interpolated to the children
				      * of the present cell and set on the
				      * terminal cells.
				      *
				      * In principle, it works as follows: if
				      * the cell pointed to by this object is
				      * terminal, then the dof values are set
				      * in the global data vector by calling
				      * the #set_dof_values# function;
				      * otherwise, the values are prolonged
				      * to each of the children and this
				      * function is called for each of them.
				      *
				      * Using the #get_interpolated_dof_values#
				      * and this function, you can compute the
				      * interpolation of a finite element
				      * function to a coarser grid by first
				      * getting the interpolated solution on a
				      * cell of the coarse grid and afterwards
				      * redistributing it using this function.
				      *
				      * Note that for continuous finite
				      * elements, calling this function affects
				      * the dof values on neighboring cells as
				      * well. It may also violate continuity
				      * requirements for hanging nodes, if
				      * neighboring cells are less refined than
				      * the present one, or if their children are
				      * less refined than the children of this
				      * cell. These requirements
				      * are not taken care of and must be
				      * enforced by the user afterwards.
				      *
				      * It is assumed that both vectors already
				      * have the right size beforehand. This
				      * function relies on the existence of a
				      * natural interpolation property of
				      * finite element spaces of a cell to
				      * its children, denoted by the
				      * prolongation matrices of finite element
				      * classes. For some elements, the spaces
				      * on coarse and fine grids are not nested,
				      * in which case the interpolation to a
				      * child is not the identity; refer to the
				      * documentation of the respective finite
				      * element class for a description of what
				      * the prolongation matrices represent in
				      * this case.
				      */
    template <typename number>
    void set_dof_values_by_interpolation (const Vector<number> &local_values,
					  Vector<number>       &values) const;
    
    				     /**
				      *  Exception
				      */
    DeclException0 (ExcNotUsefulForThisDimension);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotActive);
};





// if in optimized mode: include more templates
#ifndef DEBUG
#  include "dof_accessor.templates.h"
#endif





/*----------------------------   dof_iterator.h     ---------------------------*/
/* end of #ifndef __dof_iterator_H */
#endif
/*----------------------------   dof_iterator.h     ---------------------------*/
