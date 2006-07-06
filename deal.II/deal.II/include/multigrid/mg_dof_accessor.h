//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_dof_accessor_h
#define __deal2__mg_dof_accessor_h


#include <base/config.h>
#include <dofs/dof_accessor.h>

template <int dim> class MGDoFHandler;

template <int celldim, int dim> class MGDoFObjectAccessor;
template <int dim>              class MGDoFObjectAccessor<0, dim>;
template <int dim>              class MGDoFObjectAccessor<1, dim>;
template <int dim>              class MGDoFObjectAccessor<2, dim>;
template <int dim>              class MGDoFObjectAccessor<3, dim>;

/*!@addtogroup mg */
/*@{*/

/**
 * This is a switch class which only declares a @p typedef. It is meant to
 * determine which class aMGDoFAccessor class is to be derived
 * from. By default, <tt>MGDoFAccessor<celldim,dim></tt> derives from
 * the @p typedef in the general <tt>MGDoFObjectAccessor_Inheritance<celldim,dim></tt>
 * class, which is <tt>DoFObjectAccessor<celldim,dim></tt>,
 * but if <tt>celldim==dim</tt>, then the specialization <tt>MGDoFObjectAccessor_Inheritance<dim,dim></tt>
 * is used which declares its local type to be <tt>DoFCellAccessor<dim></tt>. Therefore,
 * the inheritance is automatically chosen to be from @p DoFCellAccessor if the
 * object under consideration has full dimension, i.e. constitutes a cell.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int celldim, int dim>
class MGDoFObjectAccessor_Inheritance 
{
  public:
				     /**
				      * Declaration of the @p typedef.
				      * See the full documentation for
				      * more information.
				      */
    typedef DoFObjectAccessor<celldim,DoFHandler<dim> > BaseClass;
};


/**
 * This is a switch class which only declares a @p typwdef. It is
 * meant to determine which class a DoFAccessor class is to be
 * derived from. By default, DoFAccessor<structdim,dim>
 * derives from the @p typedef in the general
 * <tt>DoFObjectAccessor_Inheritance<celldim,dim></tt> class, which is
 * <tt>TriaObjectAccessor<celldim,dim></tt>, but if
 * <tt>celldim==dim</tt>, then the specialization
 * <tt>TriaObjectAccessor<dim,dim></tt> is used which declares its
 * local type to be <tt>CellAccessor<dim></tt>. Therefore, the
 * inheritance is automatically chosen to be from @p CellAccessor if
 * the object under consideration has full dimension, i.e. constitutes
 * a cell.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class MGDoFObjectAccessor_Inheritance<dim,dim>
{
  public:
				     /**
				      * Declaration of the @p typedef.
				      * See the full documentation for
				      * more information.
				      */
    typedef DoFCellAccessor<DoFHandler<dim> > BaseClass;
};


/* -------------------------------------------------------------------------- */



/**
 * Define the basis for accessors to the degrees of freedom for
 * a multigrid DoF handler object.
 *
 * This class is very similar to the DoFAccessor class, and the same
 * holds for the classes derived from the present one. In essence, the
 * classes here simply extend the classes of the DoFAccessor hierarchy
 * by providing access to the multilevel degrees of freedom managed by
 * the MGDoFHandler class. See the DoFAccessor class for more
 * information on the structure of inheritance and other aspects.
 *
 * @ingroup mg
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006
 */
template <int structdim, int dim>
class MGDoFAccessor : public MGDoFObjectAccessor_Inheritance<structdim, dim>::BaseClass
{
  public:
				     /**
				      * Declare a typedef to the base
				      * class to make accessing some
				      * of the exception classes
				      * simpler.
				      */    
    typedef
    typename MGDoFObjectAccessor_Inheritance<structdim, dim>::BaseClass
    BaseClass;

				     /**
				      * A typedef necessary for the
				      * iterator classes.
				      */
    typedef MGDoFHandler<dim> AccessorData;
    
				     /**
				      * Default constructor. Creates
				      * an object that is not
				      * usable.
				      */
    MGDoFAccessor ();

				     /**
				      * Constructor.
				      */
    MGDoFAccessor (const Triangulation<dim> *tria,
		   const int                 level,
		   const int                 index,
		   const AccessorData       *local_data);

				     /**
				      * Reset the DoF handler pointer.
				      */
    void set_mg_dof_handler (MGDoFHandler<dim> *dh);

				     /**
				      * Copy operator.
				      */
    MGDoFAccessor & operator = (const MGDoFAccessor &da);
    
				     /**
				      * Return the index of the @p ith
				      * degree on the @p vertexth
				      * vertex for the level this
				      * object lives on.
				      */
    unsigned int mg_vertex_dof_index (const int level,
				      const unsigned int vertex,
				      const unsigned int i) const;

				     /**
				      * Set the index of the @p ith degree
				      * on the @p vertexth vertex to @p index
				      * for the level this object lives on.
				      */
    void set_mg_vertex_dof_index (const int level,
				  const unsigned int vertex,
				  const unsigned int i,
				  const unsigned int index) const;

				     /**
				      * Return the index of the @p ith
				      * degree of freedom of this line
				      * on the level this line lives
				      * on.
				      */
    unsigned int mg_dof_index (const int level,
			       const unsigned int i) const;

    				     /**
				      * Set the index of the @p ith degree
				      * of freedom of this line on the
				      * level this line lives on to @p index.
				      */
    void set_mg_dof_index (const int level,
			   const unsigned int i,
			   const unsigned int index) const;

				     /**
				      * Return an iterator to the @p
				      * i-th child.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<structdim, dim> >
    child (const unsigned int) const;
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const MGDoFAccessor &a);

    				     /**
				      * Exception for child classes
				      */
    DeclException0 (ExcInvalidObject);

  protected:
				     /**
				      * Store the address of the @p MGDoFHandler object
				      * to be accessed.
				      */
    MGDoFHandler<dim> *mg_dof_handler;  
};


/* -------------------------------------------------------------------------- */

/**
 * Closure class. Unused, but necessary for some template games.
 */
template<int dim>
class MGDoFObjectAccessor<0, dim>
{
  public:
    typedef void AccessorData;

				     /**
				      * Constructor. Should never be called
				      * and thus throws an error.
				      */
    MGDoFObjectAccessor (const Triangulation<dim> *,
			 const int,
			 const int,
			 const AccessorData *);
};


/**
 * Grant access to those aspects of multilevel degrees of freedom located on
 * lines that are dimension specific. See the MGDoFAccessor class for
 * more information.
 *
 * @ingroup mg
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006
 */
template <int dim>
class MGDoFObjectAccessor<1, dim> :  public MGDoFAccessor<1,dim>
{
  public:
				     /**
				      * Declare a typedef to the base
				      * class to make accessing some
				      * of the exception classes
				      * simpler.
				      */    
    typedef
    typename MGDoFObjectAccessor_Inheritance<1, dim>::BaseClass
    BaseClass;

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
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to an MGDoFHandler object.
				      */
    MGDoFObjectAccessor (const Triangulation<dim> *tria,
			 const int                 level,
			 const int                 index,
			 const AccessorData       *local_data);

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
    void get_mg_dof_indices (const int level,
			     std::vector<unsigned int> &dof_indices) const;

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
    template <typename number>
    void get_mg_dof_values (const int level,
			    const Vector<number> &values,
			    Vector<number>       &dof_values) const;
};


/**
 * Grant access to those aspects of multilevel degrees of freedom located on
 * quads that are dimension specific. See the MGDoFAccessor class for
 * more information.
 *
 * @ingroup mg
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006
 */
template <int dim>
class MGDoFObjectAccessor<2, dim> :  public MGDoFAccessor<2,dim>
{
  public:
				     /**
				      * Declare a typedef to the base
				      * class to make accessing some
				      * of the exception classes
				      * simpler.
				      */    
    typedef
    typename MGDoFObjectAccessor_Inheritance<2, dim>::BaseClass
    BaseClass;

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
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to an MGDoFHandler object.
				      */
    MGDoFObjectAccessor (const Triangulation<dim> *tria,
			 const int                 level,
			 const int                 index,
			 const AccessorData       *local_data);

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
    void get_mg_dof_indices (const int level,
			     std::vector<unsigned int> &dof_indices) const;
 
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
    template <typename number>
    void get_mg_dof_values (const int level,
			    const Vector<number> &values,
			    Vector<number>       &dof_values) const;

   				     /**
				      * Return a pointer to the @p ith line
				      * bounding this @p Quad.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
    line (const unsigned int i) const;
};


/**
 * Grant access to those aspects of multilevel degrees of freedom located on
 * hexes that are dimension specific. See the MGDoFAccessor class for
 * more information.
 *
 * @ingroup mg
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006
 */
template <int dim>
class MGDoFObjectAccessor<3, dim> :  public MGDoFAccessor<3,dim>
{
  public:
				     /**
				      * Declare a typedef to the base
				      * class to make accessing some
				      * of the exception classes
				      * simpler.
				      */    
    typedef
    typename MGDoFObjectAccessor_Inheritance<3, dim>::BaseClass
    BaseClass;

    
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
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to an MGDoFHandler object.
				      */
    MGDoFObjectAccessor (const Triangulation<dim> *tria,
			 const int                 level,
			 const int                 index,
			 const AccessorData       *local_data);

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
    void get_mg_dof_indices (const int level,
			     std::vector<unsigned int> &dof_indices) const;

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
    template <typename number>
    void get_mg_dof_values (const int level,
			    const Vector<number> &values,
			    Vector<number>       &dof_values) const;

    				     /**
				      * Return a pointer to the @p ith line
				      * bounding this @p Hex.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
    line (const unsigned int i) const;

				     /**
				      * Return a pointer to the @p ith quad
				      * bounding this @p Hex.
				      */
    TriaIterator<dim,MGDoFObjectAccessor<2, dim> >
    quad (const unsigned int i) const;
};


/**
 * Grant access to the degrees of freedom on cells, as far as this
 * isn't already covered by the MGDoFAccessor and MGDoFObjectAccessor
 * classes from which the present class is derived. In particular,
 * this function overloads functions from CellAccessor and
 * DoFCellAccessor that return iterators to other objects, such as the
 * face() or neighbor() function. Since the functions in CellAccessor
 * and DoFCellAccessor return iterators into Triangulations and
 * DoFHandlers only, we need to reimplement the functions in this
 * class to make sure we get iterators into MGDoFHandlers instead.
 *
 * @ingroup mg
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006
 */
template <int dim>
class MGDoFCellAccessor :  public MGDoFObjectAccessor<dim, dim>
{
  public:
				     /**
				      * Type of faces.
				      */
    typedef
    TriaIterator<dim, MGDoFObjectAccessor<dim-1, dim> >
    face_iterator;
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef typename MGDoFObjectAccessor<dim, dim>::AccessorData  AccessorData;
    
    				     /**
				      * Constructor
				      */
    MGDoFCellAccessor (const Triangulation<dim> *tria,
		       const int                 level,
		       const int                 index,
		       const AccessorData       *local_data);

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
    void get_mg_dof_indices (std::vector<unsigned int> &dof_indices) const;

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
    template <typename number>
    void get_mg_dof_values (const Vector<number> &values,
			    Vector<number>       &dof_values) const;

				     /**
				      * Return the @p ith neighbor as a MGDoF cell
				      * iterator. This function is needed since
				      * the neighbor function of the base
				      * class returns a cell accessor without
				      * access to the MGDoF data.
				      */
    TriaIterator<dim,MGDoFCellAccessor<dim> >
    neighbor (const unsigned int) const;

    				     /**
				      * Return the @p ith child as a MGDoF cell
				      * iterator. This function is needed since
				      * the child function of the base
				      * class returns a cell accessor without
				      * access to the DoF data.
				      */
    TriaIterator<dim,MGDoFCellAccessor<dim> >
    child (const unsigned int) const;

    				     /**
				      * Return an iterator to the @p ith face
				      * of this cell.
				      *
				      * This function is not implemented in 1D,
				      * and maps to MGDoFObjectAccessor<2, dim>::line in 2D.
				      */
    face_iterator
    face (const unsigned int i) const;

				     /**
				      * Return the result of the
				      * @p neighbor_child_on_subface
				      * function of the base class,
				      * but convert it so that one can
				      * also access the MGDoF data (the
				      * function in the base class
				      * only returns an iterator with
				      * access to the triangulation
				      * data).
				      */
    TriaIterator<dim,MGDoFCellAccessor<dim> >
    neighbor_child_on_subface (const unsigned int face_no,
                               const unsigned int subface_no) const;
};

/*@}*/


template<> MGDoFCellAccessor<1>::face_iterator
MGDoFCellAccessor<1>::face (const unsigned int i) const;
template<> MGDoFCellAccessor<2>::face_iterator
MGDoFCellAccessor<2>::face (const unsigned int i) const;
template<> MGDoFCellAccessor<3>::face_iterator
MGDoFCellAccessor<3>::face (const unsigned int i) const;


#endif
