/*----------------------------   mg_dof_accessor.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __mg_dof_accessor_H
#define __mg_dof_accessor_H
/*----------------------------   mg_dof_accessor.h     ---------------------------*/


#include <grid/dof_accessor.h>




/**
 * Define the basis for accessors to the degrees of freedom for
 * a multigrid DoF handler object.
 *
 * Note that it is allowed to construct an object of which the
 * #mg_dof_handler# pointer is a Null pointer. Such an object would
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
class MGDoFAccessor {
  public:
				     /**
				      * Constructor
				      */
    MGDoFAccessor () : mg_dof_handler(0) {
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
    MGDoFAccessor (const MGDoFHandler<dim> *mg_dof_handler) :
		    mg_dof_handler(const_cast<MGDoFHandler<dim>*>(mg_dof_handler)) {};

				     /**
				      * Reset the DoF handler pointer.
				      */
    void set_mg_dof_handler (MGDoFHandler<dim> *dh) {
      Assert (dh != 0, DoFAccessor<dim>::ExcInvalidObject());
      mg_dof_handler = dh;
    };

				     /**
				      * Copy operator.
				      */
    MGDoFAccessor<dim> & operator = (const MGDoFAccessor<dim> &da) {
      set_dof_handler (da.mg_dof_handler);
      return *this;
    };

  protected:
				     /**
				      * Store the address of the #MGDoFHandler# object
				      * to be accessed.
				      */
    MGDoFHandler<dim> *mg_dof_handler;  
};


template <int celldim, int dim, typename BaseClass>
class MGDoFObjectAccessor :  public MGDoFAccessor<dim>,
			     public DoFObjectAccessor<celldim, dim>
{};


/**
 * Closure class.
 */
template<int dim, typename BaseClass>
class MGDoFObjectAccessor<0, dim, BaseClass>
{
  public:
    typedef void* AccessorData;
    MGDoFObjectAccessor (Triangulation<dim> *,
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
 * The same concept is used with #DoFObjectAccessor<2, dim># classes etc.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim, typename BaseClass>
class MGDoFObjectAccessor<1, dim, BaseClass> :  public MGDoFAccessor<dim>,
						public DoFObjectAccessor<1, dim>
{
  public:
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef MGDoFHandler<dim> AccessorData;
    
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    MGDoFObjectAccessor ();
    
    				     /**
				      * Constructor. The #local_data#
				      * argument is assumed to be a pointer
				      * to a #MGDoFHandler<dim># object.
				      */
    MGDoFObjectAccessor (Triangulation<dim> *tria,
		       const int           level,
		       const int           index,
		       const AccessorData *local_data);
    
				     /**
				      * Return the index of the #i#th degree
				      * of freedom of this line on the level
				      * this line lives on.
				      */
    int mg_dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the #i#th degree
				      * of freedom of this line on the
				      * level this line lives on to #index#.
				      */
    void set_mg_dof_index (const unsigned int i, const int index) const;

				     /**
				      * Return the index of the #i#th degree
				      * on the #vertex#th vertex for the
				      * level this line lives on.
				      */
    int mg_vertex_dof_index (const unsigned int vertex,
			     const unsigned int i) const;

				     /**
				      * Set the index of the #i#th degree
				      * on the #vertex#th vertex to #index#
				      * for the level this line lives on.
				      */
    void set_mg_vertex_dof_index (const unsigned int vertex,
				  const unsigned int i,
				  const int          index) const;

    				     /**
				      * Return the indices of the dofs of this
				      * line in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, 
				      * dofs on line.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand. The
				      * indices refer to the local numbering
				      * for the level this line lives on.
				      */
    void get_mg_dof_indices (vector<int> &dof_indices) const;

				     /**
				      * Return the #i#th child as a MGDoF line
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a line accessor without
				      * access to the MG data.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<1, dim,BaseClass> > child (const unsigned int) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const MGDoFObjectAccessor<1, dim,BaseClass> &a);
};






/**
 * Grant access to the multilevel degrees of freedom located on quads.
 *
 * @see DoFLineAccessor
 */
template <int dim, typename BaseClass>
class MGDoFObjectAccessor<2, dim, BaseClass> :  public MGDoFAccessor<dim>,
						public DoFObjectAccessor<2, dim> {
  public:
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef MGDoFHandler<dim> AccessorData;

				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    MGDoFObjectAccessor ();

    				     /**
				      * Constructor. The #local_data#
				      * argument is assumed to be a pointer
				      * to a #DoFHandler<dim># object.
				      */
    MGDoFObjectAccessor (Triangulation<dim> *tria,
		       const int           level,
		       const int           index,
		       const AccessorData *local_data);
    
				     /**
				      * Return the index of the #i#th degree
				      * of freedom of this quad on the level
				      * this quad lives on.
				      */
    int mg_dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the #i#th degree
				      * of freedom of this quad on the
				      * level this quad lives on to #index#.
				      */
    void set_mg_dof_index (const unsigned int i, const int index) const;

				     /**
				      * Return the index of the #i#th degree
				      * on the #vertex#th vertex for the level
				      * this quad lives on.
				      */
    int mg_vertex_dof_index (const unsigned int vertex,
			     const unsigned int i) const;

				     /**
				      * Set the index of the #i#th degree
				      * on the #vertex#th vertex for the
				      * level this quad lives on to #index#.
				      */
    void set_mg_vertex_dof_index (const unsigned int vertex,
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
				      * has the right size beforehand. The
				      * indices refer to the local numbering
				      * for the level this quad lives on.
				      */
    void get_mg_dof_indices (vector<int> &dof_indices) const;

    				     /**
				      * Return a pointer to the #i#th line
				      * bounding this #Quad#.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<1, dim,TriaObjectAccessor<1, dim> > >
    line (const unsigned int i) const;

				     /**
				      * Return the #i#th child as a DoF quad
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a quad accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<2, dim, BaseClass> > child (const unsigned int) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const MGDoFObjectAccessor<2, dim,BaseClass> &a);
};




/**
 * Grant access to the multilevel degrees of freedom located on hexhedra.
 *
 * @see DoFLineAccessor
 */
template <int dim, typename BaseClass>
class MGDoFObjectAccessor<3, dim, BaseClass> :  public MGDoFAccessor<dim>,
						public DoFObjectAccessor<3, dim>
{
  public:
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef MGDoFHandler<dim> AccessorData;

				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    MGDoFObjectAccessor ();

    				     /**
				      * Constructor. The #local_data#
				      * argument is assumed to be a pointer
				      * to a #DoFHandler<dim># object.
				      */
    MGDoFObjectAccessor (Triangulation<dim> *tria,
		      const int           level,
		      const int           index,
		      const AccessorData *local_data);
    
				     /**
				      * Return the index of the #i#th degree
				      * of freedom of this hex on the level
				      * this quad lives on.
				      */
    int mg_dof_index (const unsigned int i) const;

    				     /**
				      * Set the index of the #i#th degree
				      * of freedom of this hex on the
				      * level this hex lives on to #index#.
				      */
    void set_mg_dof_index (const unsigned int i, const int index) const;

				     /**
				      * Return the index of the #i#th degree
				      * on the #vertex#th vertex for the level
				      * this hex lives on.
				      */
    int mg_vertex_dof_index (const unsigned int vertex,
			     const unsigned int i) const;

				     /**
				      * Set the index of the #i#th degree
				      * on the #vertex#th vertex for the
				      * level this hex lives on to #index#.
				      */
    void set_mg_vertex_dof_index (const unsigned int vertex,
				  const unsigned int i,
				  const int          index) const;

    				     /**
				      * Return the indices of the dofs of this
				      * hex in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand. The
				      * indices refer to the local numbering
				      * for the level this hex lives on.
				      */
    void get_mg_dof_indices (vector<int> &dof_indices) const;

    				     /**
				      * Return a pointer to the #i#th line
				      * bounding this #Hex#.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<1, dim,TriaObjectAccessor<1, dim> > >
    line (const unsigned int i) const;

				     /**
				      * Return a pointer to the #i#th quad
				      * bounding this #Hex#.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<2, dim,TriaObjectAccessor<2, dim> > >
    quad (const unsigned int i) const;

				     /**
				      * Return the #i#th child as a DoF quad
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a hex accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<3, dim,BaseClass> > child (const unsigned int) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const MGDoFObjectAccessor<3, dim,BaseClass> &a);
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
class MGDoFCellAccessor :  public MGDoFObjectAccessor<dim, dim, CellAccessor<dim> > {
  public:
				     /**
				      * Type of faces.
				      */
    typedef
    TriaIterator<dim, MGDoFObjectAccessor<dim-1, dim,TriaObjectAccessor<dim-1, dim> > >
    face_iterator;
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef typename
    MGDoFObjectAccessor<dim, dim, CellAccessor<dim> >::AccessorData
    AccessorData;
    
    				     /**
				      * Constructor
				      */
    MGDoFCellAccessor (Triangulation<dim> *tria,
		       const int           level,
		       const int           index,
		       const AccessorData *local_data) :
		    MGDoFObjectAccessor<dim, dim, CellAccessor<dim> > (tria,level,index,local_data) {};

				     /**
				      * Return the #i#th neighbor as a MGDoF cell
				      * iterator. This function is needed since
				      * the neighbor function of the base
				      * class returns a cell accessor without
				      * access to the MGDoF data.
				      */
    TriaIterator<dim,MGDoFCellAccessor<dim> > neighbor (const unsigned int) const;

    				     /**
				      * Return the #i#th child as a MGDoF cell
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a cell accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,MGDoFCellAccessor<dim> > child (const unsigned int) const;

    				     /**
				      * Return an iterator to the #i#th face
				      * of this cell.
				      *
				      * This function is not implemented in 1D,
				      * and maps to MGDoFObjectAccessor<2, dim>::line in 2D.
				      */
    face_iterator
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
				      * has the right size beforehand. The
				      * indices refer to the multilevel
				      * numbering local to the present
				      * level of this cell. The vector shall
				      * therefore have the same number of
				      * entries as there are degrees of
				      * freedom on this level.
				      */
    void get_mg_dof_values (const Vector<double> &values,
			    Vector<double>       &dof_values) const;


    				     /**
				      *  Exception
				      */
    DeclException0 (ExcNotUsefulForThisDimension);
};







/*----------------------------   mg_dof_accessor.h     ---------------------------*/
/* end of #ifndef __mg_dof_accessor_H */
#endif
/*----------------------------   mg_dof_accessor.h     ---------------------------*/
