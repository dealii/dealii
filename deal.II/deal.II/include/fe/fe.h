/*----------------------------   fe.h     ---------------------------*/
/*      $Id$                 */
#ifndef __fe_H
#define __fe_H
/*----------------------------   fe.h     ---------------------------*/

#include <lac/dfmatrix.h>
#include <grid/point.h>
#include <base/exceptions.h>
#include <grid/tria.h>



// forward declarations
template <int dim> class Boundary;
template <int dim> class FiniteElement;
template <int dim> class Quadrature;


/**
  Provide a set of flags which tells the #FEValues<>::reinit# function, which
  fields are to be updated for each cell. E.g. if you do not need the
  gradients since you want to assemble the mass matrix, you can switch that
  off. By default, all flags are off, i.e. no reinitialization will be done.
 
  A variable of this type has to be passed to the constructor of the
  #FEValues# object. You can select more than one flag by concatenation
  using the #|# (bitwise #or#) operator.
  */
enum UpdateFields {
				       /**
					* Default: update nothing.
					*/
      update_default  = 0,
				       /**
					* Compute quadrature points in real
					* space (not on unit cell).
					*/
      update_q_points = 1,
				       /**
					* Transform gradients on unit cell to
					* gradients on real cell.
					*/
      update_gradients = 2,
				       /**
					* Compute jacobian matrices of the
					* transform between unit and real cell
					* in the evaluation points.
					*/
      update_jacobians = 4,
				       /**
					* Compute the JxW values (Jacobian
					* determinant at the quadrature point
					* times the weight of this point).
					*/
      update_JxW_values = 8,
				       /**
					* Compute the points on the real cell
					* on which the ansatz functions are
					* located.
					*/
      update_ansatz_points = 16
};




/**
  Represent a finite element evaluated with a specific quadrature rule.
  This class is an optimization which avoids evaluating the shape functions
  at the quadrature points each time a quadrature takes place. Rather, the
  values and gradients (and possibly higher order derivatives in future
  versions of this library) are evaluated once and for all on the unit
  cell before doing the quadrature itself. Only the Jacobian matrix of
  the transformation from the unit cell to the real cell and the integration
  points in real space are calculated each time we move on to a new cell.

  The unit cell is defined to be the tensor product of the interval $[0,1]$
  in the present number of dimensions. In part of the literature, the convention
  is used that the unit cell be the tensor product of the interval $[-1,1]$,
  which is to distinguished properly.

  Objects of this class store a multitude of different values needed to
  do the assemblage steps on real cells rather than on the unit cell. Among
  these values are the values and gradients of the shape functions at the
  quadrature points on the real and the unit cell, the location of the
  quadrature points on the real and on the unit cell, the weights of the
  quadrature points, the Jacobian matrices of the mapping from the unit to
  the real cell at the quadrature points and so on.

  The Jacobian matrix is defined to be
  $$ J_{ij} = {d\xi_i \over dx_j} $$
  where the $\xi_i$ are the coordinates on the unit cell and the $x_i$ are
  the coordinates on the real cell.
  This is the form needed to compute the gradient on the real cell from
  the gradient on the unit cell. If we want to transform the area element
  $dx dy$ from the real to the unit cell, we have to take the determinant of
  the inverse matrix, which is the reciprocal value of the determinant of the
  matrix defined above.

  The #FEValues# object keeps track of those fields which really need to
  be computed, since the computation of the gradients of the ansatz functions
  on each real cell can be quite an expensive thing if it is not needed. The
  object knows about which fields are needed by the #UpdateFields# object
  passed through the constructor. In debug mode, the accessor functions, which
  return values from the different fields, check whether the required field
  was initialized, thus avoiding use of unitialized data.
  */
template <int dim>
class FEValues {
  public:

    
    
				     /**
				      * Number of quadrature points.
				      */
    const unsigned int n_quadrature_points;

				     /**
				      * Total number of shape functions.
				      */
    const unsigned int total_dofs;
    
				     /**
				      * Constructor. Fill all arrays with the
				      * values of the shape functions of the
				      * specified finite element using the
				      * quadrature points of the given
				      * quadrature rule.
				      *
				      * This function actually only fills
				      * the fields related to the unit face,
				      * the fields related to a real face (like
				      * gradients, true quadrature points, etc.)
				      * need to be initialized using the
				      * #reinit# function.
				      */
    FEValues (const FiniteElement<dim> &,
	      const Quadrature<dim> &,
	      const UpdateFields);

				     /**
				      * Return the value of the #i#th shape
				      * function at the #j# quadrature point.
				      */
    double shape_value (const unsigned int i,
			const unsigned int j) const;

				     /**
				      * Return a pointer to the matrix holding
				      * all values of shape functions at all
				      * integration points, on the present cell.
				      * For the format of this matrix, see the
				      * documentation for the matrix itself.
				      */
    const dFMatrix & get_shape_values () const;
    
				     /**
				      * Return the gradient of the #i#th shape
				      * function at the #j# quadrature point.
				      * If you want to get the derivative in
				      * one of the coordinate directions, use
				      * the appropriate function of the #Point#
				      * class to extract one component. Since
				      * only a reference to the gradient's value
				      * is returned, there should be no major
				      * performance drawback.
				      * The function returns the gradient on the
				      * real element, not the reference element.
				      */
    const Point<dim> & shape_grad (const unsigned int i,
				   const unsigned int j) const;

				     /** 
				      * Return a pointer to the matrix holding
				      * all gradients of shape functions at all
				      * integration points, on the present cell.
				      * For the format of this matrix, see the
				      * documentation for the matrix itself.
				      */
    const vector<vector<Point<dim> > > & get_shape_grads () const;
    
				     /**
				      * Return the position of the #i#th
				      * quadrature point in real space.
				      */
    const Point<dim> & quadrature_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of
				      * quadrature points.
				      */
    const vector<Point<dim> > & get_quadrature_points () const;

				     /**
				      * Return the point in real space where
				      * the #i#th ansatz function is located
				      * (location is in the sense of where it
				      * assumes its nominal properties, e.g. at
				      * the vertex of a cell, at the center of
				      * a line, etc).
				      *
				      * This function is needed for the
				      * interpolation problem: if we want to
				      * transfer a continuous function to a
				      * finite element function by interpolation
				      * we have to take the continuous
				      * function's value at the ansatz function
				      * locations.
				      */
    const Point<dim> & ansatz_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of points
				      * denoting the location of the ansatz
				      * functions.
				      */
    const vector<Point<dim> > & get_ansatz_points () const;
    
				     /**
				      * Return the Jacobi determinant times
				      * the weight of the #i#th quadrature
				      * point.
				      */
    double JxW (const unsigned int i) const;

				     /**
				      * Return a pointer to the array holding
				      * the JxW values at the different
				      * quadrature points.
				      */
    const vector<double> & get_JxW_values () const;
    
				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the given cell
				      * and the given finite element.
				      *
				      * This function needs a boundary object
				      * passed, since this class needs to know
				      * how to handle faces which are located
				      * on the boundary of the domain. In that
				      * case, faces may be curved and the
				      * calculation of quadrature points,
				      * gradients and the like may need
				      * additional effort, depending on the
				      * mapping from the unit to the real cell
				      * (linear mappings use straight boundary
				      * segments, but higher order elements
				      * may use other ways.)
				      */
    void reinit (const Triangulation<dim>::cell_iterator &,
		 const FiniteElement<dim> &,
		 const Boundary<dim> &);

				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The index " << arg1
		    << " is out of range, it should be less than " << arg2);
				     /**
				      * Exception
				      */
    DeclException0 (ExcAccessToUninitializedField);
				     /**
				      * Exception
				      */
    DeclException0 (ExcCannotInitializeField);
    
  private:
				     /**
				      * Store the values of the shape functions
				      * at the quadrature points. Rows in this
				      * matrix denote the values of a single
				      * shape function at the different points,
				      * columns are for a single point with the
				      * different shape functions.
				      */
    dFMatrix             shape_values;

				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points.
				      * Since unfortunately the full matrix
				      * classes of DEAL are not templated,
				      * we have to store them in an
				      * archetypic style.
				      *
				      * This field is reset each time
				      * #reinit# is called and contains the
				      * gradients on the real element, rather
				      * than on the reference element.
				      */
    vector<vector<Point<dim> > >  shape_gradients;

				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points on
				      * the unit cell.
				      * This field is set up upon construction
				      * of the object and contains the gradients
				      * on the reference element.
				      */
    vector<vector<Point<dim> > >   unit_shape_gradients;
    
				     /**
				      * Store an array of the weights of the
				      * quadrature points. This array is
				      * set up upon construction.
				      */
    vector<double>       weights;

				     /**
				      * Store an array of weights times the
				      * Jacobi determinant at the quadrature
				      * points. This function is reset each time
				      * #reinit# is called. The Jacobi determinant
				      * is actually the reciprocal value of the
				      * Jacobi matrices stored in this class,
				      * see the general documentation of this
				      * class for more information.
				      */
    vector<double>       JxW_values;

				     /**
				      * Array of quadrature points. This array
				      * is set up upon calling #reinit# and
				      * contains the quadrature points on the
				      * real element, rather than on the
				      * reference element.
				      */
    vector<Point<dim> >  quadrature_points;

				     /**
				      * Array of quadrature points in the unit
				      * cell. This array is set up upon
				      * construction and contains the quadrature
				      * points on the reference element.
				      */
    vector<Point<dim> >  unit_quadrature_points;
    
				     /**
				      * Array of points denoting the off-point
				      * of the ansatz functions. In real space
				      * (no-one seems to need the off-point
				      * on the unit cell, so no function is
				      * provided for this).
				      */
    vector<Point<dim> >  ansatz_points;
    
				     /**
				      * Store the jacobi matrices at the
				      * different quadrature points. This field
				      * is set each time #reinit# is called.
				      */
    vector<dFMatrix>     jacobi_matrices;

				     /**
				      * Store which fields are to be updated by
				      * the reinit function.
				      */
    UpdateFields         update_flags;
};




/**
  Represent a finite element evaluated with a specific quadrature rule.
  This class is an optimization which avoids evaluating the shape functions
  at the quadrature points each time a quadrature takes place. Rather, the
  values and gradients (and possibly higher order derivatives in future
  versions of this library) are evaluated once and for all on the unit
  face before doing the quadrature itself. Only the Jacobian matrix of
  the transformation from the unit face to the real face and the integration
  points in real space are calculated each time we move on to a new face.

  The unit face is defined to be the tensor product of the interval $[0,1]$
  in the present number of dimensions minus one. In part of the literature,
  the convention is used that the unit cell be the tensor product of the
  interval $[-1,1]$, which is to distinguished properly.

  This class is very similar to the #FEValues# class; see there for more
  documentation.
  */
template <int dim>
class FEFaceValues {
  public:

    
    
				     /**
				      * Number of quadrature points on
				      * the face.
				      */
    const unsigned int n_quadrature_points;

				     /**
				      * Total number of shape functions
				      * on the cell adjacent to this face.
				      * This number is not the same as the
				      * number of shape functions of which
				      * the center is located on the face.
				      */
    const unsigned int total_dofs;
    
				     /**
				      * Constructor. Fill all arrays with the
				      * values of the shape functions of the
				      * specified finite element using the
				      * quadrature points of the given
				      * quadrature rule for the face, which
				      * has a dimension one less than the
				      * cell.
				      *
				      * This function actually only fills
				      * the fields related to the unit face,
				      * the fields related to a real face (like
				      * gradients, true quadrature points, etc.)
				      * need to be initialized using the
				      * #reinit# function.
				      */
    FEFaceValues (const FiniteElement<dim> &,
		  const Quadrature<dim-1> &,
		  const UpdateFields);

				     /**
				      * Return the value of the #i#th shape
				      * function at the #j# quadrature point.
				      */
    double shape_value (const unsigned int i,
			const unsigned int j) const;

				     /**
				      * Return a pointer to the matrix holding
				      * all values of shape functions at all
				      * integration points, on the present cell.
				      * For the format of this matrix, see the
				      * documentation for the matrix itself.
				      */
    const dFMatrix & get_shape_values () const;
    
				     /**
				      * Return the gradient of the #i#th shape
				      * function at the #j# quadrature point.
				      * If you want to get the derivative in
				      * one of the coordinate directions, use
				      * the appropriate function of the #Point#
				      * class to extract one component. Since
				      * only a reference to the gradient's value
				      * is returned, there should be no major
				      * performance drawback.
				      * The function returns the gradient on the
				      * real element, not the reference element.
				      */
    const Point<dim> & shape_grad (const unsigned int i,
				   const unsigned int j) const;

				     /** 
				      * Return a pointer to the matrix holding
				      * all gradients of shape functions at all
				      * integration points, on the present cell.
				      * For the format of this matrix, see the
				      * documentation for the matrix itself.
				      */
    const vector<vector<Point<dim> > > & get_shape_grads () const;
    
				     /**
				      * Return the position of the #i#th
				      * quadrature point in real space.
				      *
				      * For curved boundary cells, using
				      * biquadratic or higher mappings
				      * of the unit cell to the real cell,
				      * these points may not be on the
				      * plane submannifold on which the
				      * vertices of the face lie.
				      */
    const Point<dim> & quadrature_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of
				      * quadrature points.
				      */
    const vector<Point<dim> > & get_quadrature_points () const;

				     /**
				      * Return the point in real space where
				      * the #i#th ansatz function is located
				      * (location is in the sense of where it
				      * assumes its nominal properties, e.g. at
				      * the vertex of a cell, at the center of
				      * a line, etc).
				      *
				      * This function is needed for the
				      * interpolation problem: if we want to
				      * transfer a continuous function to a
				      * finite element function by interpolation
				      * we have to take the continuous
				      * function's value at the ansatz function
				      * locations.
				      */
    const Point<dim> & ansatz_point (const unsigned int i) const;

				     /**
				      * Return a pointer to the vector of points
				      * denoting the location of the ansatz
				      * functions.
				      */
    const vector<Point<dim> > & get_ansatz_points () const;
    
				     /**
				      * Return the Jacobi determinant times
				      * the weight of the #i#th quadrature
				      * point.
				      */
    double JxW (const unsigned int i) const;

				     /**
				      * Return a pointer to the array holding
				      * the JxW values at the different
				      * quadrature points.
				      */
    const vector<double> & get_JxW_values () const;
    
				     /**
				      * Reinitialize the gradients, Jacobi
				      * determinants, etc for the given cell
				      * and the given finite element.
				      *
				      * The constructor needs a boundary object
				      * passed, since this class needs to know
				      * how to handle faces which are located
				      * on the boundary of the domain. In that
				      * case, faces may be curved and the
				      * calculation of quadrature points,
				      * gradients and the like may need
				      * additional effort, depending on the
				      * mapping from the unit to the real cell
				      * (linear mappings use straight boundary
				      * segments, but higher order elements
				      * may use other ways.)
				      */
    void reinit (const Triangulation<dim>::face_iterator &,
		 const FiniteElement<dim> &,
		 const Boundary<dim> &);

				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The index " << arg1
		    << " is out of range, it should be less than " << arg2);
				     /**
				      * Exception
				      */
    DeclException0 (ExcAccessToUninitializedField);
				     /**
				      * Exception
				      */
    DeclException0 (ExcCannotInitializeField);
    
  private:
				     /**
				      * Store the values of the shape functions
				      * at the quadrature points. Rows in this
				      * matrix denote the values of a single
				      * shape function at the different points,
				      * columns are for a single point with the
				      * different shape functions.
				      */
    dFMatrix             shape_values;

				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points.
				      * Since unfortunately the full matrix
				      * classes of DEAL are not templated,
				      * we have to store them in an
				      * archetypic style.
				      *
				      * This field is reset each time
				      * #reinit# is called and contains the
				      * gradients on the real element, rather
				      * than on the reference element.
				      */
    vector<vector<Point<dim> > >  shape_gradients;

				     /**
				      * Store the gradients of the shape
				      * functions at the quadrature points on
				      * the unit cell.
				      * This field is set up upon construction
				      * of the object and contains the gradients
				      * on the reference element.
				      */
    vector<vector<Point<dim> > >   unit_shape_gradients;
    
				     /**
				      * Store an array of the weights of the
				      * quadrature points. This array is
				      * set up upon construction.
				      */
    vector<double>       weights;

				     /**
				      * Store an array of weights times the
				      * Jacobi determinant at the quadrature
				      * points. This function is reset each time
				      * #reinit# is called. The Jacobi determinant
				      * is actually the reciprocal value of the
				      * Jacobi matrices stored in this class,
				      * see the general documentation of this
				      * class for more information.
				      */
    vector<double>       JxW_values;

				     /**
				      * Array of quadrature points. This array
				      * is set up upon calling #reinit# and
				      * contains the quadrature points on the
				      * real element, rather than on the
				      * reference element.
				      */
    vector<Point<dim> >  quadrature_points;

				     /**
				      * Array of quadrature points in the unit
				      * cell. This array is set up upon
				      * construction and contains the quadrature
				      * points on the reference element.
				      */
    vector<Point<dim-1> > unit_quadrature_points;
    
				     /**
				      * Array of points denoting the off-point
				      * of the ansatz functions. In real space
				      * (no-one seems to need the off-point
				      * on the unit cell, so no function is
				      * provided for this).
				      */
    vector<Point<dim> >  ansatz_points;
    
				     /**
				      * Store the jacobi matrices at the
				      * different quadrature points. This field
				      * is set each time #reinit# is called.
				      */
    vector<dFMatrix>     jacobi_matrices;

				     /**
				      * Store a pointer to the object describing
				      * the boundary of the domain.
				      */
    const Boundary<dim> &boundary;

				     /**
				      * Store which fields are to be updated by
				      * the reinit function.
				      */
    UpdateFields         update_flags;
};




/**
  Base class for finite elements in arbitrary dimensions.
  */
template <int dim>
class FiniteElementBase {
  public:
				     /**
				      * Total number of degrees of freedom
				      * on a cell. This information is
				      * redundant to some fields in the
				      * derived classes but makes
				      * writing dimension independant
				      * programs easier.
				      */
    const unsigned int total_dofs;

    				     /**
				      * Default constructor. Constructs an element
				      * which is not so useful. Checking
				      * #total_dofs# is therefore a good way to
				      * check if something went wrong.
				      */
    FiniteElementBase () :
		    total_dofs(0) {};
    
				     /**
				      * Constructor. You have to set the
				      * matrices explicitely after calling
				      * this base class' constructor.
				      */
    FiniteElementBase (const unsigned int total_dofs) :
		    total_dofs(total_dofs) {};

				     /**
				      * Copy constructor.
				      */
    FiniteElementBase (const FiniteElementBase &f);

				     /**
				      * Destructor. Only declared to have a
				      * virtual destructor which the compiler
				      * wants to have.
				      */
    virtual ~FiniteElementBase () {};
    

				     /**
				      * Return the value of the #i#th shape
				      * function at the point #p#. This function
				      * should really be pure, but then we could
				      * not make copies of a finite element
				      * object even if we did not intend to use
				      * this function. Therefore, we omit the
				      * #=0# signature and implement this function
				      * by throwing an exception.
				      * #p# is a point on the reference element.
				      */
    virtual double shape_value (const unsigned int i,
			        const Point<dim>& p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at the point #p#. This function
				      * should really be pure, but then we could
				      * not make copies of a finite element
				      * object even if we did not intend to use
				      * this function. Therefore, we omit the
				      * #=0# signature and implement this function
				      * by throwing an exception.
				      * #p# is a point on the reference element,
				      */
    virtual Point<dim> shape_grad (const unsigned int i,
				   const Point<dim>& p) const;

				     /**
				      * Return a readonly reference to the
				      * matrix which describes the transfer of a
				      * child with the given number to the
				      * mother cell.
				      */
    const dFMatrix & restrict (const unsigned int child) const;

				     /**
				      * Return a readonly reference to the
				      * matrix which describes the transfer of a
				      * mother cell to the child with the given
				      * number.
				      */
    const dFMatrix & prolongate (const unsigned int child) const;

				     /**
				      * Return a readinly reference to the
				      * matrix which describes the constraints
				      * at the interface between a refined and
				      * an unrefined cell.
				      *
				      * The matrix is obviously empty in only
				      * one space dimension, since there are no
				      * constraints then.
				      */
    const dFMatrix & constraints () const;
    
				     /**
				      * Compute the Jacobian matrix and the
				      * quadrature points as well as the ansatz
				      * function locations on the real cell in
				      * real space from the given cell
				      * and the given quadrature points on the
				      * unit cell. The Jacobian matrix is to
				      * be computed at every quadrature point.
				      * This function has to be in the finite
				      * element class, since different finite
				      * elements need different transformations
				      * of the unit cell to a real cell.
				      *
				      * The computation of the three fields may
				      * share some common code, which is why we
				      * put it in one function. However, it may
				      * not always be necessary to really
				      * compute all fields, so there are
				      * bool flags which tell the function which
				      * of the fields to actually compute.
				      *
				      * Refer to the documentation of the
				      * \Ref{FEValues} class for a definition
				      * of the Jacobi matrix.
				      *
				      * It is provided for the finite element
				      * class in one space dimension, but for
				      * higher dimensions, it depends on the
				      * present fe and needs reimplementation
				      * by the user.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      */
    virtual void fill_fe_values (const Triangulation<dim>::cell_iterator &cell,
				 const vector<Point<dim> >               &unit_points,
				 vector<dFMatrix>    &jacobians,
				 const bool           compute_jacobians,
				 vector<Point<dim> > &ansatz_points,
				 const bool           compute_ansatz_points,
				 vector<Point<dim> > &q_points,
				 const bool           compute_q_points,
				 const Boundary<dim> &boundary) const;

				     /**
				      * Return the ansatz points this FE has
				      * on a face if a cell would have the
				      * given face as a side. This function is
				      * needed for higher order elements, if
				      * we want to use curved boundary
				      * approximations. For that reason, a
				      * boundary object has to be passed.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      */
    virtual void face_ansatz_points (const Triangulation<dim>::face_iterator &face,
				     const Boundary<dim>  &boundary,
				     vector<Point<dim> >  &ansatz_points) const;
    
				     /**
				      * Comparison operator. We also check for
				      * equality of the constraint matrix,
				      * which is quite an expensive operation.
				      * Do therefore
				      * use this function with care, if possible
				      * only for debugging purposes.
				      *
				      * Since this function is not that important,
				      * we avoid an implementational question
				      * about comparing arrays and do not compare
				      * the matrix arrays #restriction# and
				      * #prolongation#.
				      */
    bool operator == (const FiniteElementBase<dim> &) const;

				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidIndex,
		    int,
		    << "Invalid index " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcPureFunctionCalled);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongFieldDimension,
		    int, int,
		    << "The field has not the assumed dimension " << arg2
		    << ", but has " << arg1 << " elements.");
  protected:
				     /**
				      * Have #N=2^dim# matrices keeping the
				      * restriction constants for the transfer
				      * of the #i#th child to the mother cell.
				      * The numbering conventions for the
				      * degree of freedom indices are descibed
				      * in the derived classes.
				      * In this matrix, the row indices belong
				      * to the destination cell, i.e. the
				      * unrefined one, while the column indices
				      * are for the refined cell's degrees of
				      * freedom.
				      *
				      * Upon assembling the transfer matrix
				      * between cells using this matrix array,
				      * zero elements in the restriction
				      * matrix are discarded and will not fill
				      * up the transfer matrix.
				      */
    dFMatrix restriction[(1<<dim)];

    				     /**
				      * Have #N=2^dim# matrices keeping the
				      * prolongation constants for the transfer
				      * of the mother cell to the #i#th child.
				      * The numbering conventions for the
				      * degree of freedom indices are descibed
				      * in the derived classes.
				      * In this matrix, the row indices belong
				      * to the destination cell, i.e. the
				      * refined one, while the column indices
				      * are for the unrefined cell's degrees of
				      * freedom.
				      *
				      * Upon assembling the transfer matrix
				      * between cells using this matrix array,
				      * zero elements in the prolongation
				      * matrix are discarded and will not fill
				      * up the transfer matrix.
				      */
    dFMatrix prolongation[(1<<dim)];

    				     /**
				      * Specify the constraints which the
				      * dofs on the two sides of a cell interface
				      * underly if the line connects two
				      * cells of which one is refined once.
				      *
				      * For further details see the general
				      * description of the derived class.
				      *
				      * This field is obviously useless in one
				      * space dimension.
				      */
    dFMatrix interface_constraints;
};






/**
  Define a finite element class.
  */
template <int dim>
class FiniteElement;



/**
  Finite Element in one dimension.

  {\bf Note on extending the finite element class}

  If you want to extend this class (not by derivation, but by adding new
  elements), you should be aware that it may be copied at places where you
  don't expect that (e.g. the #DoFHandler# class keeps a copy). You must
  thus make sure that the copy operator works correctly, in special if
  pointers are involved, copying by the default copy constructor supplied
  by the compiler will result in double deletion and maybe access to data
  elements which are no more valid.

  Consequence: make sure the copy constructor is correct.

  This class should really be a pure one, with functions having the #=0#
  signature. However, some instances of this class need to be floating around
  anyhow (e.g. the #DoFHandler# class keeps a copy, only to have the values
  of #dof_per*# available), so we do not make it pure but rather implement
  those functions which should in fact be pure to throw an error.
  */
class FiniteElement<1> : public FiniteElementBase<1> {
  public:
				     /**
				      * Number of degrees of freedom on
				      * a vertex.
				      */
    const unsigned int dofs_per_vertex;

				     /** Number of degrees of freedom
				      *  on a line.
				      */
    const unsigned int dofs_per_line;


				     /**
				      * Default constructor. The base class
				      * produces an invalid element.
				      */
    FiniteElement () :
		    dofs_per_vertex(0),
		    dofs_per_line(0) {};
    
				     /**
				      * Constructor
				      */
    FiniteElement (const unsigned int dofs_per_vertex,
		   const unsigned int dofs_per_line) :
		    FiniteElementBase<1> (2*dofs_per_vertex +
					  dofs_per_line),
  		    dofs_per_vertex(dofs_per_vertex),
		    dofs_per_line  (dofs_per_line) {};

				     /**
				      * Copy constructor
				      */
    FiniteElement (const FiniteElement<1> &fe) :
		    FiniteElementBase<1> (fe),
  		    dofs_per_vertex(fe.dofs_per_vertex),
		    dofs_per_line  (fe.dofs_per_line) {};

				     /**
				      * Same pseudo-comparison operator
				      * as for the base class.
				      */
    bool operator == (const FiniteElement<1> &f) const;

				     /**
				      * Compute the Jacobian matrix and the
				      * quadrature points as well as the ansatz
				      * function locations on the real cell in
				      * real space from the given cell
				      * and the given quadrature points on the
				      * unit cell. The Jacobian matrix is to
				      * be computed at every quadrature point.
				      *
				      * Refer to the documentation of the
				      * \Ref{FEValues} class for a definition
				      * of the Jacobi matrix.
				      *
				      * For one dimensional finite elements,
				      * these transformations are usually the
				      * same, linear ones, so we provide
				      * them in the FE<1> base class. You may,
				      * however override this implementation
				      * if you would like to use finite elements
				      * of higher than first order with
				      * non-equidistant integration points, e.g.
				      * with an exponential dependence from the
				      * distance to the origin. The standard
				      * implementation distributes the dofs on
				      * the line equidistantly.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      */
    virtual void fill_fe_values (const Triangulation<1>::cell_iterator &cell,
				 const vector<Point<1> >               &unit_points,
				 vector<dFMatrix>  &jacobians,
				 const bool         compute_jacobians,
				 vector<Point<1> > &ansatz_points,
				 const bool         compute_ansatz_points,
				 vector<Point<1> > &q_points,
				 const bool         compute_q_points,
				 const Boundary<1> &boundary) const;

				     /**
				      * Return the ansatz points this FE has
				      * on a face if a cell would have the
				      * given face as a side. This function is
				      * needed for higher order elements, if
				      * we want to use curved boundary
				      * approximations.
				      *
				      * Question: is this function useful in 1D?
				      * At present it is not implemented.
				      */
    virtual void face_ansatz_points (const Triangulation<1>::face_iterator &face,
				     const Boundary<1>  &boundary,
				     vector<Point<1> >  &ansatz_points) const;
};





/**
  Finite Element in two dimensions.

  In addition to the fields already present in 1D, a constraint matrix
  is needed in case two quads meet at a common line of which one is refined
  once more than the other one. Then there are constraints referring to the
  hanging nodes on that side of the line which is refined. These constraints
  are represented by a $n\times m$-matrix #line_constraints#, where $n$ is the
  number of degrees of freedom on the refined side (those dofs on the three
  vertices plus those on the two lines), and $m$ is that of the unrefined side
  (those dofs on the two vertices plus those on the line). The matrix is thus
  a rectangular one, being higher than wide.

  The mapping of the dofs onto the indices of the matrix is as follows:
  let $d_v$ be the number of dofs on a vertex, $d_l$ that on a line, then
  $m=0...d_v-1$ refers to the dofs on vertex zero of the unrefined line,
  $m=d_v...2d_v-1$ to those on vertex one,
  $m=2d_v...2d_v+d_l-1$ to those on the line.

  Similarly, $n=0...d_v-1$ refers to the dofs on the middle vertex
  (vertex one of child line zero, vertex zero of child line one),
  $n=d_v...d_v+d_l-1$ refers to the dofs on child line zero,
  $n=d_v+d_l...d_v+2d_l-1$ refers to the dofs on child line one.
  Please note that we do not need to reserve space for the dofs on the
  end vertices of the refined lines, since these must be mapped one-to-one
  to the appropriate dofs of the vertices of the unrefined line.

  It should be noted that it is not possible to distribute a constrained
  degree of freedom to other degrees of freedom which are themselves
  constrained. Only one level of indirection is allowed. It is not known
  at the time of this writing whether this is a constraint itself.

  If you want to extend this class (not by derivation, but by adding new
  elements), see \Ref{FiniteElement<1>}

  This class should really be a pure one, with functions having the #=0#
  signature. However, some instances of this class need to be floating around
  anyhow (e.g. the #DoFHandler# class keeps a copy, only to have the values
  of #dof_per*# available), so we do not make it pure but rather implement
  those functions which should in fact be pure to throw an error.
  */
class FiniteElement<2> : public FiniteElementBase<2> {
  public:
				     /**
				      * Number of degrees of freedom on
				      * a vertex.
				      */
    const unsigned int dofs_per_vertex;

				     /** Number of degrees of freedom
				      *  on a line.
				      */
    const unsigned int dofs_per_line;

				     /** Number of degrees of freedom
				      *  on a quad.
				      */
    const unsigned int dofs_per_quad;

				     /**
				      * Default constructor. The base class
				      * produces an invalid element.
				      */
    FiniteElement () :
		    dofs_per_vertex(0),
		    dofs_per_line(0),
		    dofs_per_quad(0) {};

				     /**
				      * Constructor
				      */
    FiniteElement (const unsigned int dofs_per_vertex,
		   const unsigned int dofs_per_line,
		   const unsigned int dofs_per_quad) :
		    FiniteElementBase<2> (4*dofs_per_vertex +
					  4*dofs_per_line +
					  dofs_per_quad),
		    dofs_per_vertex(dofs_per_vertex),
		    dofs_per_line  (dofs_per_line),
		    dofs_per_quad  (dofs_per_quad)  {};

				     /**
				      * Copy constructor
				      */
    FiniteElement (const FiniteElement<2> &fe) :
		    FiniteElementBase<2> (fe),
  		    dofs_per_vertex(fe.dofs_per_vertex),
		    dofs_per_line  (fe.dofs_per_line),
		    dofs_per_quad  (fe.dofs_per_quad) {};

				     /**
				      * Same pseudo-comparison operator
				      * as for the base class.
				      */
    bool operator == (const FiniteElement<2> &f) const;

				     /**
				      * Compute the Jacobian matrix and the
				      * quadrature points as well as the ansatz
				      * function locations on the real cell in
				      * real space from the given cell
				      * and the given quadrature points on the
				      * unit cell. The Jacobian matrix is to
				      * be computed at every quadrature point.
				      *
				      * Refer to the documentation of the
				      * \Ref{FEValues} class for a definition
				      * of the Jacobi matrix.
				      *
				      * For two dimensional finite elements,
				      * these transformations are usually
				      * dependent on the actual finite element,
				      * which is expressed by the names
				      * sub- and isoparametric elements. This
				      * function is therefore not implemented
				      * by the FE<2> base class, but is made
				      * pure virtual.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      */
    virtual void fill_fe_values (const Triangulation<2>::cell_iterator &cell,
				 const vector<Point<2> >               &unit_points,
				 vector<dFMatrix>  &jacobians,
				 const bool         compute_jacobians,
				 vector<Point<2> > &ansatz_points,
				 const bool         compute_ansatz_points,
				 vector<Point<2> > &q_points,
				 const bool         compute_q_points,
				 const Boundary<2> &boundary) const;

				     /**
				      * Return the ansatz points this FE has
				      * on a face if a cell would have the
				      * given face as a side. This function is
				      * needed for higher order elements, if
				      * we want to use curved boundary
				      * approximations.
				      *
				      * The function assumes that the fields
				      * already have the right number of
				      * elements.
				      */
    virtual void face_ansatz_points (const Triangulation<2>::face_iterator &face,
				     const Boundary<2>  &boundary,
				     vector<Point<2> >  &ansatz_points) const;
};





/*------------------------ Inline functions -----------------------------------*/



template <int dim>
inline
const dFMatrix & FEValues<dim>::get_shape_values () const {
  return shape_values;
};




template <int dim>
inline
const vector<vector<Point<dim> > > &
FEValues<dim>::get_shape_grads () const {
  Assert (update_flags | update_gradients, ExcAccessToUninitializedField());
  return shape_gradients;
};



template <int dim>
inline
const vector<Point<dim> > &
FEValues<dim>::get_quadrature_points () const {
  Assert (update_flags | update_q_points, ExcAccessToUninitializedField());
  return quadrature_points;
};



template <int dim>
inline
const vector<Point<dim> > &
FEValues<dim>::get_ansatz_points () const {
  Assert (update_flags | update_ansatz_points, ExcAccessToUninitializedField());
  return ansatz_points;
};



template <int dim>
inline
const vector<double> &
FEValues<dim>::get_JxW_values () const {
  Assert (update_flags | update_JxW_values, ExcAccessToUninitializedField());
  return JxW_values;
};



  
/*----------------------------   fe.h     ---------------------------*/
/* end of #ifndef __fe_H */
#endif
/*----------------------------   fe.h     ---------------------------*/
