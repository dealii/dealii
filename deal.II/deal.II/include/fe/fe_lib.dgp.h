//----------------------------  fe_lib.dgp.h  ---------------------------
//    Version: $Name$
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_lib.dgp.h  ---------------------------
#ifndef __deal2__fe_lib_dgp_h
#define __deal2__fe_lib_dgp_h


/*----------------------------   fe_lib.dgp.h     ---------------------------*/


#include <fe/q1_mapping.h>


/**
 * Discontinuous P1-element on hypercubes.
 *
 * This is the implementation of a linear (sic) polynomial space on a
 * d-dimensional hypercube. The shape functions are the first @p{d+1}
 * of @p{1,x,y,z}. Later on, these should be exchanged for mutually
 * orthogonal, preferably by changing the unit cell to $[-1,1]^d$.
 *
 * @author Guido Kanschat, 2000
 */
template <int dim>
class FEDG_P1 : public FEQ1Mapping<dim>
{
  public:
				     /**
				      * Constructor
				      */
    FEDG_P1 ();

				     /**
				      * Return the value of the @p{i}th shape
				      * function at point @p{p} on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the @p{i}th shape
				      * function at point @p{p} on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the @p{i}th shape function at
				      * point @p{p} on the unit cell.
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
    virtual void get_unit_support_points (typename std::vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const typename DoFHandler<dim>::cell_iterator &cell,
				     typename std::vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const typename DoFHandler<dim>::face_iterator &face,
					  typename std::vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Please note that as allowed in the
				      * documentation of the base class,
				      * this function does not implement
				      * the setting up of the local mass
				      * matrix in three space dimensions
				      * because of too high computational
				      * costs. The specified exception
				      * is thrown instead.
				      */
    virtual void get_local_mass_matrix (const typename DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * @p{restriction} and @p{prolongation}
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};




/**
 * Discontinuous P2-element on hypercubes.
 *
 * This is the implementation of a quadratic (sic) polynomial space on a
 * d-dimensional hypercube. The shape functions are those of
 * @p{1,x,y,z, x*x, x*y, x*z, y*y, y*z, z*z} applying to the space
 * dimension. Later on, these should be exchanged for mutually
 * orthogonal, preferably by changing the unit cell to $[-1,1]^d$.
 *
 * @author Guido Kanschat, 2000
 */
template <int dim>
class FEDG_P2 : public FEQ1Mapping<dim>
{
  public:
				     /**
				      * Constructor
				      */
    FEDG_P2 ();

				     /**
				      * Return the value of the @p{i}th shape
				      * function at point @p{p} on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the @p{i}th shape
				      * function at point @p{p} on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the @p{i}th shape function at
				      * point @p{p} on the unit cell.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_support_points (typename std::vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const typename DoFHandler<dim>::cell_iterator &cell,
				     typename std::vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const typename DoFHandler<dim>::face_iterator &face,
					  typename std::vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Please note that as allowed in the
				      * documentation of the base class,
				      * this function does not implement
				      * the setting up of the local mass
				      * matrix in three space dimensions
				      * because of too high computational
				      * costs. The specified exception
				      * is thrown instead.
				      */
    virtual void get_local_mass_matrix (const typename DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * @p{restriction} and @p{prolongation}
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};




/**
 * Discontinuous P3-element on hypercubes.
 *
 * This is the implementation of a cubic (sic) polynomial space on a
 * d-dimensional hypercube. The shape functions are the basis
 * polynomials spanning the space of cubic polynomials. Later on,
 * they should be exchanged for mutually orthogonal, preferably by
 * changing the unit cell to $[-1,1]^d$.
 *
 * @author Guido Kanschat, 2000
 */
template <int dim>
class FEDG_P3 : public FEQ1Mapping<dim>
{
  public:
				     /**
				      * Constructor
				      */
    FEDG_P3 ();

				     /**
				      * Return the value of the @p{i}th shape
				      * function at point @p{p} on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the @p{i}th shape
				      * function at point @p{p} on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the @p{i}th shape function at
				      * point @p{p} on the unit cell.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_support_points (typename std::vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const typename DoFHandler<dim>::cell_iterator &cell,
				     typename std::vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const typename DoFHandler<dim>::face_iterator &face,
					  typename std::vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Please note that as allowed in the
				      * documentation of the base class,
				      * this function does not implement
				      * the setting up of the local mass
				      * matrix in three space dimensions
				      * because of too high computational
				      * costs. The specified exception
				      * is thrown instead.
				      */
    virtual void get_local_mass_matrix (const typename DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * @p{restriction} and @p{prolongation}
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};




/**
 * Discontinuous P3-element on hypercubes.
 *
 * This is the implementation of a quartic (sic) polynomial space on a
 * d-dimensional hypercube. The shape functions are the basis
 * polynomials spanning the space of cubic polynomials. Later on,
 * they should be exchanged for mutually orthogonal, preferably by
 * changing the unit cell to $[-1,1]^d$.
 *
 * @author Guido Kanschat, 2000
 */
template <int dim>
class FEDG_P4 : public FEQ1Mapping<dim>
{
  public:
				     /**
				      * Constructor
				      */
    FEDG_P4 ();

				     /**
				      * Return the value of the @p{i}th shape
				      * function at point @p{p} on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the @p{i}th shape
				      * function at point @p{p} on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the @p{i}th shape function at
				      * point @p{p} on the unit cell.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_support_points (typename std::vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const typename DoFHandler<dim>::cell_iterator &cell,
				     typename std::vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const typename DoFHandler<dim>::face_iterator &face,
					  typename std::vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on what this function does.
				      *
				      * Please note that as allowed in the
				      * documentation of the base class,
				      * this function does not implement
				      * the setting up of the local mass
				      * matrix in three space dimensions
				      * because of too high computational
				      * costs. The specified exception
				      * is thrown instead.
				      */
    virtual void get_local_mass_matrix (const typename DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * @p{restriction} and @p{prolongation}
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};


/*----------------------------   fe_lib.dgp.h     ---------------------------*/

#endif
/*----------------------------   fe_lib.dgp.h     ---------------------------*/


