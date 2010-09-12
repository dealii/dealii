#ifndef __deal2__polynomials_nedelec_h
#define __deal2__polynomials_nedelec_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/tensor.h>
#include <base/point.h>
#include <base/polynomial.h>
#include <base/polynomial_space.h>
#include <base/tensor_product_polynomials.h>
#include <base/table.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN
/**
 * @addtogroup Polynomials
 * @{
 */

/**
 * This class implements the <i>H<sup>curl</sup></i>-conforming,
 * vector-valued Nédélec polynomials as proposed in the book of
 * P. Solin, K. Segeth and I. Dolezel.
 *
 * The Nédélec polynomials are constructed such that the curl
 * is in the tensor product polynomial space <i>Q<sub>k</sub></i>.
 * Therefore, the polynomial order of each component must be one
 * order higher in the corresponding two directions,
 * yielding the polynomial spaces <i>(Q<sub>k,k+1</sub>,
 * Q<sub>k+1,k</sub>)</i> and <i>(Q<sub>k,k+1,k+1</sub>,
 * Q<sub>k+1,k,k+1</sub>, Q<sub>k+1,k+1,k</sub>)</i> in 2D and 3D, resp.
 *
 * @author Markus Bürg, 2009
 */
template <int dim>
class PolynomialsNedelec
{
  public:
				     /**
				      * Constructor. Creates all basis
				      * functions for Nédélec polynomials
				      * of given degree.
				      *
				      * @arg k: the degree of the
				      * Nédélec space, which is the degree
				      * of the largest tensor product
				      * polynomial space
				      * <i>Q<sub>k</sub></i> contained.
				      */
    PolynomialsNedelec (const unsigned int k);
    
				     /**
				      * Computes the value and the
				      * first and second derivatives
				      * of each Nédélec
				      * polynomial at @p unit_point.
				      *
				      * The size of the vectors must
				      * either be zero or equal
				      * <tt>n()</tt>.  In the
				      * first case, the function will
				      * not compute these values.
				      *
				      * If you need values or
				      * derivatives of all tensor
				      * product polynomials then use
				      * this function, rather than
				      * using any of the
				      * <tt>compute_value</tt>,
				      * <tt>compute_grad</tt> or
				      * <tt>compute_grad_grad</tt>
				      * functions, see below, in a
				      * loop over all tensor product
				      * polynomials.
				      */
    void compute (const Point<dim> &unit_point, std::vector<Tensor<1, dim> > &values, std::vector<Tensor<2, dim> > &grads, std::vector<Tensor<3, dim> > &grad_grads) const;

				     /**
				      * Returns the number of Nédélec
					* polynomials.
				      */
    unsigned int n () const;
    
				     /**
				      * Returns the degree of the Nédélec
				      * space, which is one less than
				      * the highest polynomial degree.
				      */
    unsigned int degree () const;
    
				     /**
				      * Return the number of
				      * polynomials in the space
				      * <TT>N(degree)</tt> without
				      * requiring to build an object
				      * of PolynomialsNedelec. This is
				      * required by the FiniteElement
				      * classes.
				      */
    static unsigned int compute_n_pols (unsigned int degree);
    
  private:
				     /**
				      * The degree of this object as
				      * given to the constructor.
				      */
    const unsigned int my_degree;
    
				     /**
				      * An object representing the
				      * polynomial space for a single
				      * component. We can re-use it by
				      * rotating the coordinates of
				      * the evaluation point.
				      */
    const AnisotropicPolynomials<dim> polynomial_space;

				     /**
				      * Number of Nédélec polynomials.
				      */
    const unsigned int n_pols;

				     /**
				      * A static member function that
				      * creates the polynomial space
				      * we use to initialize the
				      * #polynomial_space member
				      * variable.
				      */
    static std::vector<std::vector< Polynomials::Polynomial< double > > > create_polynomials (const unsigned int k);
};

/** @} */

template <int dim>
inline unsigned int PolynomialsNedelec<dim>::n () const
{
  return n_pols;
}

template <int dim>
inline unsigned int PolynomialsNedelec<dim>::degree () const
{
  return my_degree;
}
DEAL_II_NAMESPACE_CLOSE

#endif
