/*----------------------------   dof_iterator.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dof_iterator_H
#define __dof_iterator_H
/*----------------------------   dof_iterator.h     ---------------------------*/


template <int dim> class Triangulation;
template <int dim> class DoFHandler;
class dVector;


#include <grid/tria_accessor.h>
#include <vector>



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
				      */
    DoFAccessor (DoFHandler<dim> *dof_handler) :
		    dof_handler(dof_handler) {};

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
				      * Exception for child classes
				      */
    DeclException3 (ExcInvalidIndex,
		    int,
		    int,
		    int,
		    << "Invalid index " << arg1
		    << ", index must be between " << arg2
		    << " and " << arg3 << ".");
				     /**
				      * Exception
				      */
    DeclException0 (ExcVectorNotEmpty);
				     /**
				      * Exception
				      */
    DeclException0 (ExcVectorDoesNotMatch);

  protected:
				     /**
				      * Store the address of the #DoFHandler# object
				      * to be accessed.
				      */
    DoFHandler<dim> *dof_handler;  
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
 * #DoFLineAccessor<1># and #CellAccessor<1># (the #DoFLineAccessor<dim># cannot
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
 * is a #LineAccessor<dim>#, i.e. am accessor to lines without the missing
 * functionality needed for cells (neighbors, etc.).
 *
 * This way we can declare a #DoFCellAccessor# in one dimension by deriving
 * from #DoFLineAccessor<1,CellAccessor<1> >#, thus getting the cell
 * functionality through the #DoFLineAccessor# instead of through a virtual
 * multiple inheritance of #DoFLineAccessor# and #CellAccessor<1>#.
 *
 * The same concept is used with #DoFQuadAccessor# classes etc.
 */
template <int dim, typename BaseClass>
class DoFLineAccessor :  public DoFAccessor<dim>, public BaseClass {
  public:
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFLineAccessor ();
    
    				     /**
				      * Constructor. The #local_data#
				      * argument is assumed to be a pointer
				      * to a #DoFHandler<dim># object.
				      */
    DoFLineAccessor (Triangulation<dim> *tria,
		     const int           level,
		     const int           index,
		     const void         *local_data) :
		    DoFAccessor<dim> ((DoFHandler<dim>*)local_data),
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
    TriaIterator<dim,DoFLineAccessor<dim,BaseClass> > child (const unsigned int) const;

				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFLineAccessor<dim,BaseClass> &a);
};






/**
 * Grant access to the degrees of freedom located on quads.
 *
 * @see DoFLineAccessor
 */
template <int dim, typename BaseClass>
class DoFQuadAccessor :  public DoFAccessor<dim>, public BaseClass {
  public:
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFQuadAccessor ();

    				     /**
				      * Constructor. The #local_data#
				      * argument is assumed to be a pointer
				      * to a #DoFHandler<dim># object.
				      */
    DoFQuadAccessor (Triangulation<dim> *tria,
		     const int           level,
		     const int           index,
		     const void         *local_data) :
		    DoFAccessor<dim> ((DoFHandler<dim>*)local_data),
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
    TriaIterator<dim,DoFLineAccessor<dim,LineAccessor<dim> > >
    line (const unsigned int i) const;

				     /**
				      * Return the #i#th child as a DoF quad
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a quad accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,DoFQuadAccessor<dim,BaseClass> > child (const unsigned int) const;

				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFQuadAccessor<dim,BaseClass> &a);
};








/**
 * Intermediate, "typedef"-class, not for public use.
 *
 * Rationale for the declaration of members for this class: gcc 2.8 has a bug
 * when deriving from explicitely specialized classes which materializes in
 * the calculation of wrong addresses of member variables. By declaring the
 * general template of #DoFSubstructAccessor# to have the same object layout as
 * the specialized versions (using the same base classes), we fool the compiler,
 * which still looks in the wrong place for the addresses but finds the
 * right information. This way, at least ot works.
 *
 * Insert a guard, however, in the constructor to avoid that anyone (including
 * the compiler) happens to use this class.
 */
template <int dim>
class DoFSubstructAccessor : public DoFAccessor<dim>,
			     public TriaAccessor<dim> {
  public:
    DoFSubstructAccessor () {
      Assert (false, ExcInternalError());
    };

    DeclException0 (ExcInternalError);
};




/**
 * Intermediate, "typedef"-class, not for public use.
 *
 * \subsection{Rationale}
 *
 * This class is only a wrapper class used to do kind of a typedef
 * with template parameters. This class and #DoFSubstructAccessor<2>#
 * wrap the following names:
 * \begin{verbatim}
 *   DoFSubstructAccessor<1> := DoFLineAccessor<1,CellAccessor<1> >;
 *   DoFSubstructAccessor<2> := DoFQuadAccessor<2,CellAccessor<2> >;
 * \end{verbatim}
 * We do this rather complex (and needless, provided C++ the needed constructs!)
 * class hierarchy manipulation, since this way we can declare and implement
 * the \Ref{DoFCellAccessor} dimension independent as an inheritance from
 * #DoFSubstructAccessor<dim>#. If we had not declared these
 * types, we would have to write two class declarations, one for
 * #DoFCellAccessor<1>#, derived from #DoFLineAccessor<1,CellAccessor<1> >#
 * and one for #DoFCellAccessor<2>#, derived from
 * #DoFQuadAccessor<2,CellAccessor<2> >#.
 */
template <>
class DoFSubstructAccessor<1> :  public DoFLineAccessor<1,CellAccessor<1> > {
  public:
    				     /**
				      * Constructor
				      */
    DoFSubstructAccessor (Triangulation<1> *tria,
			  const int         level,
			  const int         index,
			  const void       *local_data) :
		    DoFLineAccessor<1,CellAccessor<1> > (tria,level,index,local_data) {};
				     // do this here, since this way the
				     // CellAccessor has the possibility to know
				     // what a face_iterator is. Otherwise
				     // it would have to ask the DoFHandler<dim>
				     // which would need another big include
				     // file and maybe cyclic interdependence
    typedef void * face_iterator;
};



/**
 * Intermediate, "typedef"-class, not for public use.
 *
 * @see DoFSubstructAccessor<1>
 */
template <>
class DoFSubstructAccessor<2> : public DoFQuadAccessor<2,CellAccessor<2> > {
  public:
    				     /**
				      * Constructor
				      */
    DoFSubstructAccessor (Triangulation<2> *tria,
			  const int         level,
			  const int         index,
			  const void       *local_data) :
		    DoFQuadAccessor<2,CellAccessor<2> > (tria,level,index,local_data) {};
				     // do this here, since this way the
				     // CellAccessor has the possibility to know
				     // what a face_iterator is. Otherwise
				     // it would have to ask the DoFHandler<dim>
				     // which would need another big include
				     // file and maybe cyclic interdependence
    typedef TriaIterator<2,DoFLineAccessor<2,LineAccessor<2> > > face_iterator;
};






/**
 * Grant access to the degrees of freedom on a cell. In fact, since all
 * access to the degrees of freedom has been enabled by the classes
 * #DoFLineAccessor<1># and #DoFQuadAccessor<2># for the space dimension
 * one and two, respectively, this class only collects the pieces
 * together by deriving from the appropriate #DoF*Accessor# and the
 * right #CellAccessor<dim># and finally adding two functions which give
 * access to the neighbors and children as #DoFCellAccessor# objects
 * rather than #CellAccessor# objects (the latter function was inherited
 * from the #CellAccessor<dim># class).
 */
template <int dim>
class DoFCellAccessor :  public DoFSubstructAccessor<dim> {
  public:
    				     /**
				      * Constructor
				      */
    DoFCellAccessor (Triangulation<dim> *tria,
		     const int           level,
		     const int           index,
		     const void         *local_data) :
		     DoFSubstructAccessor<dim> (tria,level,index,local_data) {};

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
				      * and maps to DoFQuadAccessor::line in 2D.
				      */
    typename DoFSubstructAccessor<dim>::face_iterator
    face (const unsigned int i) const;

    				     /**
				      * Return the value of the given vector
				      * restricted to the dofs of this
				      * cell in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand.
				      */
    void get_dof_values (const dVector  &values,
			 vector<double> &dof_values) const;

    				     /**
				      *  Exception
				      */
    DeclException0 (ExcNotUsefulForThisDimension);
};










/*----------------------------   dof_iterator.h     ---------------------------*/
/* end of #ifndef __dof_iterator_H */
#endif
/*----------------------------   dof_iterator.h     ---------------------------*/
