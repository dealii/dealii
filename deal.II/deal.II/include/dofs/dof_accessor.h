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
#ifndef __deal2__dof_accessor_h
#define __deal2__dof_accessor_h


#include <base/config.h>
#include <grid/tria_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/hp_dof_handler.h>
#include <vector>


template <typename number> class FullMatrix;
template <typename number> class SparseMatrix;
template <typename number> class Vector;

//TODO: Find a way to avoid the forward declarations and get rid of dof_handler include files.

template <int celldim, class DH> class DoFObjectAccessor;

template <int dim, typename Accessor> class TriaRawIterator;



// note: the file dof_accessor.templates.h is included at the end of
// this file.  this includes a lot of templates and thus makes
// compilation slower, but at the same time allows for more aggressive
// inlining and thus faster code.


/**
 * This is a switch class which only declares a @p typedef. It is meant to
 * determine which class a DoFObjectAccessor class is to be derived
 * from. By default, <tt>DoFObjectAccessor<celldim,dim></tt> derives from
 * the @p typedef in the general <tt>DoFObjectAccessor_Inheritance<celldim,dim></tt>
 * class, which is <tt>TriaObjectAccessor<celldim,dim></tt>,
 * but if <tt>celldim==dim</tt>, then the specialization <tt>DoFObjectAccessor_Inheritance<dim,dim></tt>
 * is used which declares its local type to be <tt>CellAccessor<dim></tt>. Therefore,
 * the inheritance is automatically chosen to be from CellAccessor if the
 * object under consideration has full dimension, i.e. constitutes a cell.
 *
 * @ingroup dofs
 * @ingroup Accessors 
 * @author Wolfgang Bangerth, 1999
 */
template <int celldim, int dim>
class DoFObjectAccessor_Inheritance 
{
  public:
				     /**
				      * Declaration of the @p typedef.
				      * See the full documentation for
				      * more information.
				      */
    typedef TriaObjectAccessor<celldim,dim> BaseClass;
};


/**
 * This is a switch class which only declares a @p typedef. It is meant to
 * determine which class a DoFObjectAccessor class is to be derived
 * from. By default, <tt>DoFObjectAccessor<celldim,dim></tt> derives from
 * the @p typedef in the general <tt>DoFObjectAccessor_Inheritance<celldim,dim></tt>
 * class, which is <tt>TriaObjectAccessor<celldim,dim></tt>,
 * but if <tt>celldim==dim</tt>, then the specialization <tt>DoFObjectAccessor_Inheritance<dim,dim></tt>
 * is used which declares its local type to be <tt>CellAccessor<dim></tt>. Therefore,
 * the inheritance is automatically chosen to be from CellAccessor if the
 * object under consideration has full dimension, i.e. constitutes a cell.
 *
 * @ingroup dofs
 * @ingroup Accessors 
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class DoFObjectAccessor_Inheritance<dim,dim>
{
  public:
				     /**
				      * Declaration of the @p typedef.
				      * See the full documentation for
				      * more information.
				      */
    typedef CellAccessor<dim> BaseClass;
};


/* -------------------------------------------------------------------------- */



/**
 * Define the base class for accessors to the degrees of
 * freedom. Accessors are used to, well, access the data that pertains
 * to edges, faces, and cells of a triangulation. The concept is
 * explained in more detail in connection to @ref Iterators.
 *
 * This class follows mainly the route laid out by the accessor library
 * declared in the triangulation library (TriaAccessor). It enables
 * the user to access the degrees of freedom on the lines (there are similar
 * versions for the DoFs on quads, etc), where the dimension of the underlying
 * triangulation does not really matter (i.e. this accessor works with the
 * lines in 1D-, 2D-, etc dimensions).
 *
 * The first template argument of this class refers to the structural
 * dimension of the thing accessed: it is 1 for lines, 2 for quads, or
 * 3 for hexes. The second argument denotes the type of DoF handler we
 * should work on. It can either be ::DoFHandler<dim> or
 * hp::DoFHandler<dim>. The space dimension of the object we are to
 * work on is automatically extracted from the DH template argument.
 *
 * Depending on whether the structural dimension of the object
 * accessed equals the space dimension on which the DoF handler object
 * operates, this class is derived from CellAccessor or
 * TriaObjectAccessor. This means that, for example accessors to quads
 * in 2d have access to all the mesh aspects of cells, whereas
 * accessors to quads in 3d can only access things that make sense for
 * faces.
 *
 *
 * <h3>Usage</h3>
 *
 * Usage is best to happen through the typedefs to the various kinds
 * of iterators provided by the DoFHandler and hp::DoFHandler classes,
 * since they are more secure to changes in the class naming and
 * template interface as well as providing easier typing (much less
 * complicated names!).
 * 
 *
 * <h3>Inheritance</h3>
 *
 * If the structural dimension given by the first template argument
 * equals the dimension of the DoFHandler (given as the second
 * template argument), then we are obviously dealing with cells,
 * rather than lower-dimensional objects. In that case, inheritance is
 * from CellAccessor, to provide access to all the cell specific
 * information afforded by that class. Otherwise, i.e. for
 * lower-dimensional objects, inheritance is from TriaObjectAccessor.
 *
 * There are classes DoFObjectAccessor<1>, DoFObjectAccessor<2>, and
 * DoFObjectAccessor<3> derived from the present class that provide
 * things that are specific to objects of a particular
 * dimensionality. For example, DoFObjectAccessor<2> allows access to
 * the DoF iterators pointing to the lines that bound that
 * quadrilateral, which is obviously something that a line can't
 * do. Consequently, these few operations are declared in these
 * derived classes.
 *
 * In addition, there is a DoFCellAccessor class that provides the
 * equivalent to the CellAccessor class.
 *
 * @ingroup dofs
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006
 */
template <int structdim, class DH>
class DoFAccessor : public DoFObjectAccessor_Inheritance<structdim, DH::dimension>::BaseClass
{
  public:
				     /**
				      * Declare a typedef to the base
				      * class to make accessing some
				      * of the exception classes
				      * simpler.
				      */
    typedef
    typename DoFObjectAccessor_Inheritance<structdim, DH::dimension>::BaseClass
    BaseClass;

				     /**
				      * Data type passed by the iterator class.
				      */
    typedef DH AccessorData;

				     /**
				      * Default constructor. Provides
				      * an accessor that can't be
				      * used.
				      */
    DoFAccessor ();
    
				     /**
				      * Constructor
				      */
    DoFAccessor (const Triangulation<DH::dimension> *tria,
		 const int                 level,
		 const int                 index,
		 const DH                 *local_data);

				     /**
				      * Reset the DoF handler pointer.
				      */
    void set_dof_handler (DH *dh);

				     /**
				      * Return a handle on the
				      * DoFHandler object which we
				      * are using.
				      */
    const DH &
    get_dof_handler () const;

				     /**
				      * Copy operator.
				      */
    DoFAccessor<structdim,DH> &
    operator = (const DoFAccessor<structdim,DH> &da);
    
				     /**
				      * Implement the copy operator needed
				      * for the iterator classes.
				      */
    void copy_from (const DoFAccessor<structdim, DH> &a);

				     /**
				      * Return an iterator pointing to
				      * the the @p c-th child.
				      */
    TriaIterator<DH::dimension,DoFObjectAccessor<structdim,DH> >
    child (const unsigned int c) const;

				     /**
				      * Index of the <i>i</i> degree
				      * on the @p vertexth vertex.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    unsigned int vertex_dof_index (const unsigned int vertex,
				   const unsigned int i,
				   const unsigned int fe_index = DH::default_fe_index) const;

				     /**
				      * Set the index of the <i>i</i> degree
				      * on the @p vertex-th vertex to @p index.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    void set_vertex_dof_index (const unsigned int vertex,
			       const unsigned int i,
			       const unsigned int index,
			       const unsigned int fe_index = DH::default_fe_index) const;

				     /**
				      * Index of the <i>i</i>th degree
				      * of freedom of this object.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    unsigned int dof_index (const unsigned int i,
			    const unsigned int fe_index = DH::default_fe_index) const;

    				     /**
				      * Set the index of the
				      * <i>i</i>th degree of freedom
				      * of this object to @p index.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      */
    void set_dof_index (const unsigned int i,
			const unsigned int index,
			const unsigned int fe_index = DH::default_fe_index) const;

                                     /**
                                      * Return the number of finite
                                      * elements that are active on a
                                      * given object.
                                      *
                                      * For non-hp DoFHandler objects,
                                      * the answer is of course always
                                      * one. However, for
                                      * hp::DoFHandler objects, this
                                      * isn't the case: If this is a
                                      * cell, the answer is of course
                                      * one. If it is a face, the
                                      * answer may be one or two,
                                      * depending on whether the two
                                      * adjacent cells use the same
                                      * finite element or not. If it
                                      * is an edge in 3d, the possible
                                      * return value may be one or any
                                      * other value larger than that.
                                      */
    unsigned int
    n_active_fe_indices () const;

				     /**
				      * Return the @p n-th active fe
				      * index on this object. For
				      * cells and all non-hp objects,
				      * there is only a single active
				      * fe index, so the argument must
				      * be equal to zero. For
				      * lower-dimensional hp objects,
				      * there are
				      * n_active_fe_indices() active
				      * finite elements, and this
				      * function can be queried for
				      * their indices.
				      */
    unsigned int
    nth_active_fe_index (const unsigned int n) const;
    
				     /**
				      * Return true if the finite
				      * element with given index is
				      * active on the present
				      * object. For non-hp DoF
				      * accessors, this is of course
				      * the case only if @p fe_index
				      * equals zero. For cells, it is
				      * the case if @p fe_index equals
				      * active_fe_index() of this
				      * cell. For faces and other
				      * lower-dimensional objects,
				      * there may be more than one @p
				      * fe_index that are active on
				      * any given object (see
				      * n_active_fe_indices()).
				      */
    bool
    fe_index_is_active (const unsigned int fe_index) const;

                                     /**
				      * Exceptions for child classes
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcInvalidObject);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcVectorNotEmpty);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcVectorDoesNotMatch);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcMatrixDoesNotMatch);
				     /**
				      * A function has been called for
				      * a cell which should be active,
				      * but is refined. @ref GlossActive
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcNotActive);
				     /**
				      * Exception
				      * 
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcCantCompareIterators);

  protected:    

				     /**
				      * Store the address of the DoFHandler object
				      * to be accessed.
				      */
    DH *dof_handler;

				     /**
				      *  Compare for equality.            
				      */
    bool operator == (const DoFAccessor &) const;
	
				     /**
				      * Compare for inequality.
				      */
    bool operator != (const DoFAccessor &) const;

                                     /**
                                      * Iterator classes need to be friends
                                      * because they need to access operator==
                                      * and operator!=.
                                      */
    template <int, typename> friend class TriaRawIterator;
};


/* -------------------------------------------------------------------------- */


/**
 * Closure class. Unused, but necessary for some template constructs.
 * @ingroup dofs
 * @ingroup Accessors
 */
template <class DH>
class DoFObjectAccessor<0, DH> : public DoFAccessor<0, DH>
{
  public:
    				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;

    typedef void AccessorData;

				     /**
				      * Constructor. Should never be called
				      * and thus throws an error.
				      */
    DoFObjectAccessor (Triangulation<dim> *,
		       const int,
		       const int,
		       const AccessorData *)
      {
	Assert (false, ExcInternalError());
      }
};



/**
 * Grant access to those aspects of degrees of freedom located on
 * quads that are dimension specific. See the DoFAccessor class for
 * more information.
 *
 * @ingroup dofs
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006
 */
template <class DH>
class DoFObjectAccessor<1, DH> : public DoFAccessor<1,DH>
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;
    
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef DH AccessorData;

				     /**
				      * Declare base class as a local typedef
				      * for simpler access.
				      */
    typedef DoFAccessor<1,DH> BaseClass;

				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFObjectAccessor ();
    
    				     /**
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to a DoFHandler object.
				      */
    DoFObjectAccessor (const Triangulation<DH::dimension> *tria,
		       const int                 level,
		       const int                 index,
		       const AccessorData       *local_data);    

    				     /**
				      * Return the indices of the dofs of this
				      * object in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * The vector has to have the
				      * right size before being passed
				      * to this function.
				      *
				      * This function is most often
				      * used on active objects (edges,
				      * faces, cells). It can be used
				      * on non-active objects as well
				      * (i.e. objects that have
				      * children), but only if the
				      * finite element under
				      * consideration has degrees of
				      * freedom exclusively on
				      * vertices. Otherwise, the
				      * function doesn't make much
				      * sense, since for example
				      * inactive edges do not have
				      * degrees of freedom associated
				      * with them at all.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      *
				      * For cells, there is only a
				      * single possible finite element
				      * index (namely the one for that
				      * cell, returned by
				      * <code>cell-@>active_fe_index</code>. Consequently,
				      * the derived DoFCellAccessor
				      * class has an overloaded
				      * version of this function that
				      * calls the present function
				      * with
				      * <code>cell-@>active_fe_index</code>
				      * as last argument.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices,
			  const unsigned int fe_index = DH::default_fe_index) const;
};


/**
 * Grant access to those aspects of degrees of freedom located on
 * quads that are dimension specific. See the DoFAccessor class for
 * more information.
 *
 * @ingroup dofs
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006
 */
template <class DH>
class DoFObjectAccessor<2, DH> : public DoFAccessor<2,DH>
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;

				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef DH AccessorData;

				     /**
				      * Declare base class as a local typedef
				      * for simpler access.
				      */
    typedef DoFAccessor<2,DH> BaseClass;
    
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFObjectAccessor ();

    				     /**
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to a DoFHandler object.
				      */
    DoFObjectAccessor (const Triangulation<DH::dimension> *tria,
		       const int                 level,
		       const int                 index,
		       const AccessorData       *local_data);
    
    				     /**
				      * Return the indices of the dofs of this
				      * object in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * The vector has to have the
				      * right size before being passed
				      * to this function.
				      *
				      * This function is most often
				      * used on active objects (edges,
				      * faces, cells). It can be used
				      * on non-active objects as well
				      * (i.e. objects that have
				      * children), but only if the
				      * finite element under
				      * consideration has degrees of
				      * freedom exclusively on
				      * vertices. Otherwise, the
				      * function doesn't make much
				      * sense, since for example
				      * inactive edges do not have
				      * degrees of freedom associated
				      * with them at all.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      *
				      * For cells, there is only a
				      * single possible finite element
				      * index (namely the one for that
				      * cell, returned by
				      * <code>cell-@>active_fe_index</code>. Consequently,
				      * the derived DoFCellAccessor
				      * class has an overloaded
				      * version of this function that
				      * calls the present function
				      * with
				      * <code>cell-@>active_fe_index</code>
				      * as last argument.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices,
			  const unsigned int fe_index = DH::default_fe_index) const;

    				     /**
				      *  Return a pointer to the @p ith line
				      *  bounding this @p Quad.
				      */
    TriaIterator<DH::dimension,DoFObjectAccessor<1, DH> >
    line (const unsigned int i) const;
};


/**
 * Grant access to those aspects of degrees of freedom located on
 * hexes that are dimension specific. See the DoFAccessor class for
 * more information.
 *
 * @ingroup dofs
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998, 2006
 */
template <class DH>
class DoFObjectAccessor<3, DH> : public DoFAccessor<3,DH>
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;
    
				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef DH AccessorData;

				     /**
				      * Declare base class as a local typedef
				      * for simpler access.
				      */
    typedef DoFAccessor<1,DH> BaseClass;
    
				     /**
				      * Default constructor, unused thus
				      * not implemented.
				      */
    DoFObjectAccessor ();

    				     /**
				      * Constructor. The @p local_data
				      * argument is assumed to be a pointer
				      * to a DoFHandler object.
				      */
    DoFObjectAccessor (const Triangulation<DH::dimension> *tria,
		       const int                 level,
		       const int                 index,
		       const AccessorData       *local_data);
    
    				     /**
				      * Return the indices of the dofs of this
				      * object in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * The vector has to have the
				      * right size before being passed
				      * to this function.
				      *
				      * This function is most often
				      * used on active objects (edges,
				      * faces, cells). It can be used
				      * on non-active objects as well
				      * (i.e. objects that have
				      * children), but only if the
				      * finite element under
				      * consideration has degrees of
				      * freedom exclusively on
				      * vertices. Otherwise, the
				      * function doesn't make much
				      * sense, since for example
				      * inactive edges do not have
				      * degrees of freedom associated
				      * with them at all.
				      *
				      * The last argument denotes the
				      * finite element index. For the
				      * standard ::DoFHandler class,
				      * this value must be equal to
				      * its default value since that
				      * class only supports the same
				      * finite element on all cells
				      * anyway.
				      *
				      * However, for hp objects
				      * (i.e. the hp::DoFHandler
				      * class), different finite
				      * element objects may be used on
				      * different cells. On faces
				      * between two cells, as well as
				      * vertices, there may therefore
				      * be two sets of degrees of
				      * freedom, one for each of the
				      * finite elements used on the
				      * adjacent cells. In order to
				      * specify which set of degrees
				      * of freedom to work on, the
				      * last argument is used to
				      * disambiguate. Finally, if this
				      * function is called for a cell
				      * object, there can only be a
				      * single set of degrees of
				      * freedom, and fe_index has to
				      * match the result of
				      * active_fe_index().
				      *
				      * For cells, there is only a
				      * single possible finite element
				      * index (namely the one for that
				      * cell, returned by
				      * <code>cell-@>active_fe_index</code>. Consequently,
				      * the derived DoFCellAccessor
				      * class has an overloaded
				      * version of this function that
				      * calls the present function
				      * with
				      * <code>cell-@>active_fe_index</code>
				      * as last argument.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices,
			  const unsigned int fe_index = DH::default_fe_index) const;

    				     /**
				      *  Return a pointer to the @p ith line
				      *  bounding this @p Hex.
				      */
    TriaIterator<DH::dimension,DoFObjectAccessor<1, DH> >
    line (const unsigned int i) const;

    				     /**
				      *  Return a pointer to the @p ith quad
				      *  bounding this @p Hex.
				      */
    TriaIterator<DH::dimension,DoFObjectAccessor<2, DH> >
    quad (const unsigned int i) const;
};



/**
 * Grant access to the degrees of freedom on a cell. In fact, since all
 * access to the degrees of freedom has been enabled by the classes
 * <tt>DoFObjectAccessor<1, 1></tt> and <tt>DoFObjectAccessor<2, 2></tt> for the space dimension
 * one and two, respectively, this class only collects the pieces
 * together by deriving from the appropriate <tt>DoF*Accessor</tt> and the
 * right <tt>CellAccessor<dim></tt> and finally adding two functions which give
 * access to the neighbors and children as DoFCellAccessor objects
 * rather than CellAccessor objects (the latter function was inherited
 * from the <tt>CellAccessor<dim></tt> class).
 *
 * Note that since for the class we derive from, i.e. <tt>DoFObjectAccessor<dim,dim></tt>,
 * the two template parameters are equal, the base class is actually derived from
 * CellAccessor, which makes the functions of this class available to the
 * DoFCellAccessor class as well.
 *
 * @ingroup dofs
 * @ingroup Accessors
 * @author Wolfgang Bangerth, 1998
 */
template <class DH>
class DoFCellAccessor :  public DoFObjectAccessor<DH::dimension,DH>
{
  public:
				     /**
				      * Extract dimension from DH.
				      */
    static const unsigned int dim = DH::dimension;
    
				     /**
				      * Declare the data type that
				      * this accessor class expects to
				      * get passed from the iterator
				      * classes.
				      */
    typedef typename DoFObjectAccessor<DH::dimension,DH>::AccessorData AccessorData;

				     /**
				      * Declare a typedef to the base
				      * class to make accessing some
				      * of the exception classes
				      * simpler.
				      */
    typedef DoFObjectAccessor<DH::dimension,DH> BaseClass;
    
    				     /**
				      * Constructor
				      */
    DoFCellAccessor (const Triangulation<DH::dimension> *tria,
		     const int                 level,
		     const int                 index,
		     const AccessorData       *local_data);

				     /**
				      * Return the @p ith neighbor as
				      * a DoF cell iterator. This
				      * function is needed since the
				      * neighbor function of the base
				      * class returns a cell accessor
				      * without access to the DoF
				      * data.
				      */
    TriaIterator<DH::dimension,DoFCellAccessor<DH> >
    neighbor (const unsigned int) const;

    				     /**
				      * Return the @p ith child as a
				      * DoF cell iterator. This
				      * function is needed since the
				      * child function of the base
				      * class returns a cell accessor
				      * without access to the DoF
				      * data.
				      */
    TriaIterator<DH::dimension,DoFCellAccessor<DH> >
    child (const unsigned int) const;

    				     /**
				      * Return an iterator to the @p ith face
				      * of this cell.
				      *
				      * This function is not implemented in
				      * 1D, and maps to DoFObjectAccessor<2,
				      * dim>::line in 2D.
				      */
    TriaIterator<DH::dimension, DoFObjectAccessor<dim-1, DH> >
    face (const unsigned int i) const;

                                     /**
				      * Return the result of the
				      * @p neighbor_child_on_subface
				      * function of the base class,
				      * but convert it so that one can
				      * also access the DoF data (the
				      * function in the base class
				      * only returns an iterator with
				      * access to the triangulation
				      * data).
				      */
    TriaIterator<DH::dimension,DoFCellAccessor<DH> >
    neighbor_child_on_subface (const unsigned int face_no,
                               const unsigned int subface_no) const;

    				     /**
				      * Return the values of the given vector
				      * restricted to the dofs of this
				      * cell in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * The vector has to have the
				      * right size before being passed
				      * to this function. This
				      * function is only callable for
				      * active cells.
				      *
				      * The input vector may be either a
				      * <tt>Vector<float></tt>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class InputVector, typename number>
    void get_dof_values (const InputVector &values,
			 Vector<number>    &local_values) const;

				     /**
				      * This function is the counterpart to
				      * get_dof_values(): it takes a vector
				      * of values for the degrees of freedom
				      * of the cell pointed to by this iterator
				      * and writes these values into the global
				      * data vector @p values. This function
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
				      * The vector has to have the
				      * right size before being passed
				      * to this function.
				      *
				      * The output vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class OutputVector, typename number>
    void set_dof_values (const Vector<number> &local_values,
			 OutputVector         &values) const;
    
                                     /**
				      * Return the interpolation of
				      * the given finite element
				      * function to the present
				      * cell. In the simplest case,
				      * the cell is a terminal one,
				      * i.e. has no children; then,
				      * the returned value is the
				      * vector of nodal values on that
				      * cell. You could then as well
				      * get the desired values through
				      * the @p get_dof_values
				      * function In the other case,
				      * when the cell has children, we
				      * use the restriction matrices
				      * provided by the finite element
				      * class to compute the
				      * interpolation from the
				      * children to the present cell.
				      *
				      * It is assumed that both
				      * vectors already have the right
				      * size beforehand. This function
				      * assumes the existence of an
				      * interpolation from child cells
				      * to the mother cell, denoted by
				      * the restriction matrices of
				      * the finite element class.
				      * Futhermore, this interpolation
				      * should be possible for each
				      * child alone, i.e.  it should
				      * be possible to compute the
				      * restriction by writing values
			              * obtained from each child
				      * directly into the output
				      * vector, without, for example,
				      * computing an average over all
				      * children.  These properties,
				      * however, do not exist for all
				      * elements; an example is the
				      * DG(0) element, which
				      * represents piecewise constant
				      * elements: for these the
				      * restriction to mother cell
				      * could be the average of the
				      * values of the children, maybe
				      * weighted by the measure of
				      * each child. It is not yet
				      * decided what the this function
				      * does in these cases.
				      *
				      * Unlike the get_dof_values()
				      * function, this function is
				      * associated to cells rather
				      * than to lines, quads, and
				      * hexes, since interpolation is
				      * presently only provided for
				      * cells by the finite element
				      * objects.
				      */
    template <class InputVector, typename number>
    void get_interpolated_dof_values (const InputVector &values,
				      Vector<number>    &interpolated_values) const;

				     /**
				      * This, again, is the
				      * counterpart to
				      * get_interpolated_dof_values():
				      * you specify the dof values on
				      * a cell and these are
				      * interpolated to the children
				      * of the present cell and set on
				      * the terminal cells.
				      *
				      * In principle, it works as
				      * follows: if the cell pointed
				      * to by this object is terminal,
				      * then the dof values are set in
				      * the global data vector by
				      * calling the set_dof_values()
				      * function; otherwise, the
				      * values are prolonged to each
				      * of the children and this
				      * function is called for each of
				      * them.
				      *
				      * Using the
				      * get_interpolated_dof_values()
				      * and this function, you can
				      * compute the interpolation of a
				      * finite element function to a
				      * coarser grid by first getting
				      * the interpolated solution on a
				      * cell of the coarse grid and
				      * afterwards redistributing it
				      * using this function.
				      *
				      * Note that for continuous
				      * finite elements, calling this
				      * function affects the dof
				      * values on neighboring cells as
				      * well. It may also violate
				      * continuity requirements for
				      * hanging nodes, if neighboring
				      * cells are less refined than
				      * the present one, or if their
				      * children are less refined than
				      * the children of this
				      * cell. These requirements are
				      * not taken care of and must be
				      * enforced by the user
				      * afterwards.
				      *
				      * It is assumed that both
				      * vectors already have the right
				      * size beforehand. This function
				      * relies on the existence of a
				      * natural interpolation property
				      * of finite element spaces of a
				      * cell to its children, denoted
				      * by the prolongation matrices
				      * of finite element classes. For
				      * some elements, the spaces on
				      * coarse and fine grids are not
				      * nested, in which case the
				      * interpolation to a child is
				      * not the identity; refer to the
				      * documentation of the
				      * respective finite element
				      * class for a description of
				      * what the prolongation matrices
				      * represent in this case.
				      *
				      * Unlike the set_dof_values()
				      * function, this function is
				      * associated to cells rather
				      * than to lines, quads, and
				      * hexes, since interpolation is
				      * presently only provided for
				      * cells by the finite element
				      * objects.
				      *
				      * The output vector may be either a
				      * Vector<float>,
				      * Vector<double>, or a
				      * BlockVector<double>, or a
				      * PETSc vector if deal.II is compiled to
				      * support PETSc. It is in the
				      * responsibility of the caller to assure
				      * that the types of the numbers stored
				      * in input and output vectors are
				      * compatible and with similar accuracy.
				      */
    template <class OutputVector, typename number>
    void set_dof_values_by_interpolation (const Vector<number> &local_values,
					  OutputVector         &values) const;

    				     /**
				      * Return the indices of the dofs of this
				      * quad in the standard ordering: dofs
				      * on vertex 0, dofs on vertex 1, etc,
				      * dofs on line 0, dofs on line 1, etc,
				      * dofs on quad 0, etc.
				      *
				      * It is assumed that the vector already
				      * has the right size beforehand.
				      *
				      * This function reimplements the
				      * same function in the base
				      * class. The functions in the
				      * base classes are available for
				      * all geometric objects,
				      * i.e. even in 3d they can be
				      * used to access the dof indices
				      * of edges, for example. On the
				      * other hand, the most common
				      * case is clearly the use on
				      * cells, which is why we cache
				      * the array for each cell, but
				      * not edge. To retrieve the
				      * cached values, rather than
				      * collect the necessary
				      * information every time, this
				      * function overwrites the one in
				      * the base class.
				      *
				      * This function is most often
				      * used on active objects (edges,
				      * faces, cells). It can be used
				      * on non-active objects as well
				      * (i.e. objects that have
				      * children), but only if the
				      * finite element under
				      * consideration has degrees of
				      * freedom exclusively on
				      * vertices. Otherwise, the
				      * function doesn't make much
				      * sense, since for example
				      * inactive edges do not have
				      * degrees of freedom associated
				      * with them at all.
				      */
    void get_dof_indices (std::vector<unsigned int> &dof_indices) const;

                                     /**
                                      * Return the finite element that
                                      * is used on the cell pointed to
                                      * by this iterator. For non-hp
                                      * DoF handlers, this is of
                                      * course always the same
                                      * element, independent of the
                                      * cell we are presently on, but
                                      * for hp DoF handlers, this may
                                      * change from cell to cell.
                                      */
    const FiniteElement<DH::dimension> &
    get_fe () const;

                                     /**
				      *  Returns the index inside the
				      *  hp::FECollection of the FiniteElement
				      *  used for this cell.
				      */
    unsigned int active_fe_index () const;

                                     /**
				      *  Sets the index of the FiniteElement used for
				      *  this cell.
				      */
    void set_active_fe_index (const unsigned int i);

				     /**
				      * Distribute a local (cell
				      * based) vector to a global one
				      * by mapping the local numbering
				      * of the degrees of freedom to
				      * the global one and entering
				      * the local values into the
				      * global vector.
				      *
				      * The elements are
				      * <em>added</em> up to the
				      * elements in the global vector,
				      * rather than just set, since
				      * this is usually what one
				      * wants.
				      */
    template <typename number, typename OutputVector>
    void
    distribute_local_to_global (const Vector<number> &local_source,
                                OutputVector         &global_destination) const;

				     /**
				      * This function does much the
				      * same as the
				      * <tt>distribute_local_to_global(Vector,Vector)</tt>
				      * function, but operates on
				      * matrices instead of
				      * vectors. If the matrix type is
				      * a sparse matrix then it is
				      * supposed to have non-zero
				      * entry slots where required.
				      */
    template <typename number, typename OutputMatrix>
    void
    distribute_local_to_global (const FullMatrix<number> &local_source,
                                OutputMatrix             &global_destination) const;
    
  private:
				     /**
				      * Update the cache in which we
				      * store the dof indices of this
				      * cell.
				      */
    void update_cell_dof_indices_cache () const;

				     /**
				      * Make the DoFHandler class a
				      * friend so that it can call the
				      * update_cell_dof_indices_cache()
				      * function
				      */
    template <int> friend class DoFHandler;
};


/* -------------- declaration of explicit specializations ------------- */

template <>
TriaIterator<1, DoFObjectAccessor<0,DoFHandler<1> > >
DoFCellAccessor<DoFHandler<1> >::face (const unsigned int) const;
template <>
TriaIterator<2, DoFObjectAccessor<1,DoFHandler<2> > >
DoFCellAccessor<DoFHandler<2> >::face (const unsigned int i) const;
template <>
TriaIterator<3, DoFObjectAccessor<2,DoFHandler<3> > >
DoFCellAccessor<DoFHandler<3> >::face (const unsigned int i) const;


template <>
const FiniteElement<1> &
DoFCellAccessor<hp::DoFHandler<1> >::get_fe () const;
template <>
const FiniteElement<2> &
DoFCellAccessor<hp::DoFHandler<2> >::get_fe () const;
template <>
const FiniteElement<3> &
DoFCellAccessor<hp::DoFHandler<3> >::get_fe () const;

template <>
unsigned int
DoFCellAccessor<hp::DoFHandler<1> >::active_fe_index () const;
template <>
unsigned int
DoFCellAccessor<hp::DoFHandler<2> >::active_fe_index () const;
template <>
unsigned int
DoFCellAccessor<hp::DoFHandler<3> >::active_fe_index () const;

template <>
void
DoFCellAccessor<hp::DoFHandler<1> >::set_active_fe_index (const unsigned int i);
template <>
void
DoFCellAccessor<hp::DoFHandler<2> >::set_active_fe_index (const unsigned int i);
template <>
void
DoFCellAccessor<hp::DoFHandler<3> >::set_active_fe_index (const unsigned int i);

template <>
void
DoFCellAccessor<DoFHandler<1> >::update_cell_dof_indices_cache () const;
template <>
void
DoFCellAccessor<DoFHandler<2> >::update_cell_dof_indices_cache () const;
template <>
void
DoFCellAccessor<DoFHandler<3> >::update_cell_dof_indices_cache () const;

template <>
void
DoFCellAccessor<hp::DoFHandler<1> >::update_cell_dof_indices_cache () const;
template <>
void
DoFCellAccessor<hp::DoFHandler<2> >::update_cell_dof_indices_cache () const;
template <>
void
DoFCellAccessor<hp::DoFHandler<3> >::update_cell_dof_indices_cache () const;



// include more templates
#include "dof_accessor.templates.h"


/*----------------------------   dof_iterator.h     ---------------------------*/

#endif
/*----------------------------   dof_iterator.h     ---------------------------*/
