/*----------------------------   fe_lib.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __fe_lib_lagrange_H
#define __fe_lib_lagrange_H
/*----------------------------   fe_lib.h     ---------------------------*/


#include <fe/fe_linear_mapping.h>



/**
 * Define a (bi-, tri-, etc)linear finite element in #dim# space dimensions,
 * along with (bi-, tri-)linear (therefore isoparametric) transforms from the
 * unit cell to the real cell.
 *
 * The linear, isoparametric mapping from a point $\vec \xi$ on the unit cell
 * to a point $\vec x$ on the real cell is defined as
 * $$ \vec x(\vec \xi)  = \sum_j {\vec p_j} N_j(\xi) $$
 * where $\vec p_j$ is the vector to the $j$th corner point of the cell in
 * real space and $N_j(\vec \xi)$ is the value of the basis function associated
 * with the $j$th corner point, on the unit cell at point $\vec \xi$. The sum
 * over $j$ runs over all corner points.
 *
 * @author Wolfgang Bangerth, 1998
 */

template <int dim>
class FELinear : public FELinearMapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FELinear ();
  protected:
				     /**
				      * Constructor that is called by the
				      * constructor of the derived
				      * #FEDGLinear# class.
				      * It uses  no dofs in the vertices and 
				      * $2^d$ dofs per cell. No constraint
				      * matrices are build.
				      * For more detail see class #FEDGLinear#.
				      */
    FELinear (const int);
  public:

				     /**
				      * Declare a static function which returns
				      * the number of degrees of freedom per
				      * vertex, line, face, etc. This function
				      * is used by the constructors, but is
				      * mainly needed for the composed finite
				      * elements, which assemble a FE object
				      * from several other FEs. See the
				      * documentation for the #FiniteElement#
				      * class for more information on this
				      * subject.
				      */
    static const FiniteElementData<dim> get_fe_data ();
    
				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the #i#th shape function at
				      * point #p# on the unit cell.
				      *
				      * For linear elements, all second
				      * derivatives on the unit cell are zero.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_support_points (vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * #restriction# and #prolongation#
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};




/**
 * Define a (bi-, tri-, etc)quadratic finite element in #dim# space dimensions.
 * A linear (subparametric) mapping from the unit cell
 * to the real cell is implemented.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FEQuadraticSub : public FELinearMapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FEQuadraticSub ();
  protected:
				     /**
				      * Constructor that is called by the
				      * constructor of the derived
				      * #FEDGQuadraticSub# class.
				      * It uses no dofs in the vertices, no
				      * dofs in the lines and 
				      * $3^d$ dofs per cell. No constraint
				      * matrices are build.
				      * For more detail see class 
				      * #FEDGQuadraticSub#.
				      */
    FEQuadraticSub (const int);
    
  public:
				     /**
				      * Declare a static function which returns
				      * the number of degrees of freedom per
				      * vertex, line, face, etc. This function
				      * is used by the constructors, but is
				      * mainly needed for the composed finite
				      * elements, which assemble a FE object
				      * from several other FEs. See the
				      * documentation for the #FiniteElement#
				      * class for more information on this
				      * subject.
				      */
    static const FiniteElementData<dim> get_fe_data ();

				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the #i#th shape function at
				      * point #p# on the unit cell.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_support_points (vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * #restriction# and #prolongation#
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};




/**
 * Define a (bi-, tri-, etc)cubic finite element in #dim# space dimensions.
 * A linear (subparametric) mapping from the unit cell
 * to the real cell is implemented.
 *
 * The numbering of degrees of freedom in one spatial dimension is as follows:
 * \begin{verbatim}
 *   0--2--3--1
 * \end{verbatim}
 *
 * The numbering of degrees of freedom in two spatial dimension is as follows:
 * \begin{verbatim}
 *   3--8--9--2
 *   |        |
 *   11 15 14 7
 *   |        |
 *   10 12 13 6
 *   |        |
 *   0--4--5--1
 * \end{verbatim}
 * Note the reverse ordering of degrees of freedom on the left and upper
 * line and the counterclockwise numbering of the interior degrees of
 * freedom.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FECubicSub : public FELinearMapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FECubicSub ();

  protected:
				     /**
				      * Constructor that is called by the
				      * constructor of the derived
				      * #FEDGCubicSub# class.
				      * It uses  no dofs in the vertices and 
				      * $4^d$ dofs per cell. No constraint
				      * matrices are build.
				      * For more detail see class
				      * #FEDGCubicSub#.
				      */
    FECubicSub (const int);

  public:

				     /**
				      * Declare a static function which returns
				      * the number of degrees of freedom per
				      * vertex, line, face, etc. This function
				      * is used by the constructors, but is
				      * mainly needed for the composed finite
				      * elements, which assemble a FE object
				      * from several other FEs. See the
				      * documentation for the #FiniteElement#
				      * class for more information on this
				      * subject.
				      */
    static const FiniteElementData<dim> get_fe_data ();

				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the #i#th shape function at
				      * point #p# on the unit cell.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_support_points (vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * #restriction# and #prolongation#
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};



/**
 * Define a (bi-, tri-, etc)quartic finite element in #dim# space dimensions.
 * A linear (subparametric) mapping from the unit cell
 * to the real cell is implemented.
 *
 * The numbering of degrees of freedom in one spatial dimension is as follows:
 * \begin{verbatim}
 *   0--2--3--4--1
 * \end{verbatim}
 *
 * The numbering of degrees of freedom in two spatial dimension is as follows:
 * \begin{verbatim}
 *   3--10-11-12-2
 *   |           |
 *   15 19 22 18 9
 *   |           |
 *   14 23 24 21 8
 *   |           |
 *   13 16 20 17 7
 *   |           |
 *   0--4--5--6--1
 * \end{verbatim}
 * Note the reverse ordering of degrees of freedom on the left and upper
 * line and the numbering of the interior degrees of
 * freedom.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FEQuarticSub : public FELinearMapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FEQuarticSub ();

  protected:
				     /**
				      * Constructor that is called by the
				      * constructor of the derived
				      * #FEDGQuarticSub# class.
				      * It uses  no dofs in the vertices and 
				      * $5^d$ dofs per cell. No constraint
				      * matrices are build.
				      * For more detail see class
				      * #FEDGQuarticSub#.
				      */
    FEQuarticSub (const int);

  public:

				     /**
				      * Declare a static function which returns
				      * the number of degrees of freedom per
				      * vertex, line, face, etc. This function
				      * is used by the constructors, but is
				      * mainly needed for the composed finite
				      * elements, which assemble a FE object
				      * from several other FEs. See the
				      * documentation for the #FiniteElement#
				      * class for more information on this
				      * subject.
				      */
    static const FiniteElementData<dim> get_fe_data ();

				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the #i#th shape function at
				      * point #p# on the unit cell.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_support_points (vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * #restriction# and #prolongation#
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};








/*----------------------------   fe_lib.h     ---------------------------*/
/* end of #ifndef __fe_lib_lagrange_H */
#endif
/*----------------------------   fe_lib.h     ---------------------------*/


