//----------------------  tensor_product_polynomials.h  -------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------  tensor_product_polynomials.h  -------------
#ifndef __deal2__tensor_product_polynomials_h
#define __deal2__tensor_product_polynomials_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/tensor.h>
#include <base/point.h>
#include <base/polynomial.h>

#include <vector>


/**
 * Tensor product of given polynomials.
 *
 * Given a vector of <i>n</i> one-dimensional polynomials
 * <i>P<sub>1</sub></i> to <i>P<sub>n</sub></i>, this class generates
 * <i>n<sup>dim</sup></i> polynomials of the form
 * <i>Q<sub>ijk</sub>(x,y,z) =
 * P<sub>i</sub>(x)P<sub>j</sub>(y)P<sub>k</sub>(z)</i>. If the base
 * polynomials are mutually orthogonal on the interval [-1,1] or
 * [0,1], then the tensor product polynomials are orthogonal on
 * [-1,1]<sup>dim</sup> or [0,1]<sup>dim</sup>, respectively.
 *
 * Indexing is as follows: the order of dim-dimensional polynomials is
 * x-coordinates running fastest, then y-coordinate, etc. The first
 * few polynomials are thus <i>P<sub>1</sub>(x)P<sub>1</sub>(y),
 * P<sub>2</sub>(x)P<sub>1</sub>(y), P<sub>3</sub>(x)P<sub>1</sub>(y),
 * ..., P<sub>1</sub>(x)P<sub>2</sub>(y),
 * P<sub>2</sub>(x)P<sub>2</sub>(y), P<sub>3</sub>(x)P<sub>2</sub>(y),
 * ...</i> and likewise in 3d.
 * 
 * @author Ralf Hartmann, Guido Kanschat, 2000, Wolfgang Bangerth 2003
 */
template <int dim>
class TensorProductPolynomials
{
  public:
				     /**
				      * Constructor. <tt>pols</tt> is
				      * a vector of objects that
				      * should be derived or otherwise
				      * convertible to one-dimensional
				      * polynomial objects. It will be
				      * copied element by element into
				      * a private variable.
				      */
    template <class Pol>
    TensorProductPolynomials (const std::vector<Pol> &pols);

				     /**
				      * Computes the value and the
				      * first and second derivatives
				      * of each tensor product
				      * polynomial at <tt>unit_point</tt>.
				      *
				      * The size of the vectors must
				      * either be equal 0 or equal
				      * n(). In the first case, the
				      * function will not compute
				      * these values.
				      *
				      * If you need values or
				      * derivatives of all tensor
				      * product polynomials then use
				      * this function, rather than
				      * using any of the
				      * compute_value(),
				      * compute_grad() or
				      * compute_grad_grad()
				      * functions, see below, in a
				      * loop over all tensor product
				      * polynomials.
				      */
    void compute (const Point<dim>            &unit_point,
                  std::vector<double>         &values,
                  std::vector<Tensor<1,dim> > &grads,
                  std::vector<Tensor<2,dim> > &grad_grads) const;
    
				     /**
				      * Computes the value of the
				      * <tt>i</tt>th tensor product
				      * polynomial at
				      * <tt>unit_point</tt>. Here <tt>i</tt> is
				      * given in tensor product
				      * numbering.
				      *
				      * Note, that using this function
				      * within a loop over all tensor
				      * product polynomials is not
				      * efficient, because then each
				      * point value of the underlying
				      * (one-dimensional) polynomials
				      * is (unnecessarily) computed
				      * several times.  Instead use
				      * the compute() function with
				      * <tt>values.size()==</tt>n()
				      * to get the point values of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    double compute_value (const unsigned int i,
			  const Point<dim> &p) const;

				     /**
				      * Computes the grad of the
				      * <tt>i</tt>th tensor product
				      * polynomial at
				      * <tt>unit_point</tt>. Here <tt>i</tt> is
				      * given in tensor product
				      * numbering.
				      *
				      * Note, that using this function
				      * within a loop over all tensor
				      * product polynomials is not
				      * efficient, because then each
				      * derivative value of the
				      * underlying (one-dimensional)
				      * polynomials is (unnecessarily)
				      * computed several times.
				      * Instead use the compute()
				      * function, see above, with
				      * <tt>grads.size()==</tt>n()
				      * to get the point value of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    Tensor<1,dim> compute_grad (const unsigned int i,
				const Point<dim> &p) const;

				     /**
				      * Computes the second
				      * derivative (grad_grad) of the
				      * <tt>i</tt>th tensor product
				      * polynomial at
				      * <tt>unit_point</tt>. Here <tt>i</tt> is
				      * given in tensor product
				      * numbering.
				      *
				      * Note, that using this function
				      * within a loop over all tensor
				      * product polynomials is not
				      * efficient, because then each
				      * derivative value of the
				      * underlying (one-dimensional)
				      * polynomials is (unnecessarily)
				      * computed several times.
				      * Instead use the compute()
				      * function, see above, with
				      * <tt>grad_grads.size()==</tt>n()
				      * to get the point value of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    Tensor<2,dim> compute_grad_grad (const unsigned int i,
                                     const Point<dim> &p) const;

				     /**
				      * Returns the number of tensor
				      * product polynomials. For <i>n</i>
				      * 1d polynomials this is <i>n<sup>dim</sup></i>.
				      */
    unsigned int n () const;

				     /**
				      * Exception.
				      */
    DeclException3 (ExcDimensionMismatch2,
		    int, int, int,
		    << "Dimension " << arg1 << " not equal to " << arg2 << " nor to " << arg3);

	    
  private:
				     /**
				      * Copy of the vector <tt>pols</tt> of
				      * polynomials given to the
				      * constructor.
				      */
    std::vector<Polynomials::Polynomial<double> > polynomials;

				     /**
				      * Number of tensor product
				      * polynomials. See n().
				      */
    unsigned int n_tensor_pols;

                                     /**
                                      * Each tensor product polynomial
                                      * <i>i</i> is a product of
                                      * one-dimensional polynomials in
                                      * each space direction. Compute
                                      * the indices of these
                                      * one-dimensional polynomials
                                      * for each space direction,
                                      * given the index <i>i</i>.
                                      */
    void compute_index (const unsigned int i,
                        unsigned int       (&indices)[dim]) const;
    
				     /**
				      * Computes
				      * <i>x<sup>dim</sup></i> for
				      * unsigned int <i>x</i>. Used in
				      * the constructor.
				      */
    static
    unsigned int x_to_the_dim (const unsigned int x);
};




/**
 * Anisotropic tensor product of given polynomials.
 *
 * Given one-dimensional polynomials <tt>Px1</tt>, <tt>Px2</tt>, ... in
 * x-direction, <tt>Py1</tt>, <tt>Py2</tt>, ... in y-direction, and so on, this
 * class generates polynomials of the form  <i>Q<sub>ijk</sub>(x,y,z) =
 * Pxi(x)Pyj(y)Pzk(z)</i>. If the base polynomials are mutually
 * orthogonal on the interval $[-1,1]$ or $[0,d]$, then the tensor
 * product polynomials are orthogonal on $[-1,1]^d$ or $[0,1]^d$,
 * respectively.
 *
 * Indexing is as follows: the order of dim-dimensional polynomials
 * is x-coordinates running fastest, then y-coordinate, etc. The first
 * few polynomials are thus <tt>Px1(x)Py1(y)</tt>, <tt>Px2(x)Py1(y)</tt>,
 * <tt>Px3(x)Py1(y)</tt>, ..., <tt>Px1(x)Py2(y)</tt>, <tt>Px2(x)Py2(y)</tt>,
 * <tt>Px3(x)Py2(y)</tt>, ..., and likewise in 3d.
 * 
 * @author Wolfgang Bangerth 2003
 */
template <int dim>
class AnisotropicPolynomials
{
  public:
				     /**
				      * Constructor. <tt>pols</tt> is a
				      * table of one-dimensional
				      * polynomials. The number of
				      * rows in this table should be
				      * equal to the space dimension,
				      * with the elements of each row
				      * giving the polynomials that
				      * shall be used in this
				      * particular coordinate
				      * direction. These polynomials
				      * may vary between coordinates,
				      * as well as their number.
				      */
    AnisotropicPolynomials (const std::vector<std::vector<Polynomials::Polynomial<double> > > &pols);

				     /**
				      * Computes the value and the
				      * first and second derivatives
				      * of each tensor product
				      * polynomial at <tt>unit_point</tt>.
				      *
				      * The size of the vectors must
				      * either be equal <tt>0</tt> or equal
				      * <tt>n_tensor_pols</tt>.  In the
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
    void compute (const Point<dim>            &unit_point,
                  std::vector<double>         &values,
                  std::vector<Tensor<1,dim> > &grads,
                  std::vector<Tensor<2,dim> > &grad_grads) const;
    
				     /**
				      * Computes the value of the
				      * <tt>i</tt>th tensor product
				      * polynomial at
				      * <tt>unit_point</tt>. Here <tt>i</tt> is
				      * given in tensor product
				      * numbering.
				      *
				      * Note, that using this function
				      * within a loop over all tensor
				      * product polynomials is not
				      * efficient, because then each
				      * point value of the underlying
				      * (one-dimensional) polynomials
				      * is (unnecessarily) computed
				      * several times.  Instead use
				      * the <tt>compute</tt> function, see
				      * above, with
				      * <tt>values.size()==n_tensor_pols</tt>
				      * to get the point values of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    double compute_value (const unsigned int i,
			  const Point<dim> &p) const;

				     /**
				      * Computes the grad of the
				      * <tt>i</tt>th tensor product
				      * polynomial at
				      * <tt>unit_point</tt>. Here <tt>i</tt> is
				      * given in tensor product
				      * numbering.
				      *
				      * Note, that using this function
				      * within a loop over all tensor
				      * product polynomials is not
				      * efficient, because then each
				      * derivative value of the
				      * underlying (one-dimensional)
				      * polynomials is (unnecessarily)
				      * computed several times.
				      * Instead use the <tt>compute</tt>
				      * function, see above, with
				      * <tt>grads.size()==n_tensor_pols</tt>
				      * to get the point value of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    Tensor<1,dim> compute_grad (const unsigned int i,
				const Point<dim> &p) const;

				     /**
				      * Computes the second
				      * derivative (grad_grad) of the
				      * <tt>i</tt>th tensor product
				      * polynomial at
				      * <tt>unit_point</tt>. Here <tt>i</tt> is
				      * given in tensor product
				      * numbering.
				      *
				      * Note, that using this function
				      * within a loop over all tensor
				      * product polynomials is not
				      * efficient, because then each
				      * derivative value of the
				      * underlying (one-dimensional)
				      * polynomials is (unnecessarily)
				      * computed several times.
				      * Instead use the <tt>compute</tt>
				      * function, see above, with
				      * <tt>grad_grads.size()==n_tensor_pols</tt>
				      * to get the point value of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    Tensor<2,dim> compute_grad_grad (const unsigned int i,
                                     const Point<dim> &p) const;

				     /**
				      * Returns the number of tensor
				      * product polynomials. It is the
				      * product of the number of
				      * polynomials in each coordinate
				      * direction.
				      */
    unsigned int n () const;

				     /**
				      * Exception.
				      */
    DeclException3 (ExcDimensionMismatch2,
		    int, int, int,
		    << "Dimension " << arg1 << " not equal to " << arg2 << " nor to " << arg3);
                                     /**
                                      * Exception
                                      */
    DeclException1 (ExcInvalidDim,
                    int,
                    << "The number of rows in this table must be equal to the "
                    << "space dimension, but is " << arg1);
	    
  private:
				     /**
				      * Copy of the vector <tt>pols</tt> of
				      * polynomials given to the
				      * constructor.
				      */
    std::vector<std::vector<Polynomials::Polynomial<double> > > polynomials;

				     /**
				      * Number of tensor product
				      * polynomials. This is
				      * <tt>Nx*Ny*Nz</tt>, or with terms
				      * dropped if the number of space
				      * dimensions is less than 3.
				      */
    unsigned int n_tensor_pols;

                                     /**
                                      * Each tensor product polynomial
                                      * @þ{i} is a product of
                                      * one-dimensional polynomials in
                                      * each space direction. Compute
                                      * the indices of these
                                      * one-dimensional polynomials
                                      * for each space direction,
                                      * given the index <tt>i</tt>.
                                      */
    void compute_index (const unsigned int i,
                        unsigned int       (&indices)[dim]) const;
    
				     /**
				      * Given the input to the
				      * constructor, compute
				      * <tt>n_tensor_pols</tt>.
				      */
    static
    unsigned int
    get_n_tensor_pols (const std::vector<std::vector<Polynomials::Polynomial<double> > > &pols);
};


/// @if NoDoc

/* -------------- declaration of explicit specializations --- */

template <>
void
TensorProductPolynomials<1>::compute_index(const unsigned int n,
                                           unsigned int      (&index)[1]) const;
template <>
void
TensorProductPolynomials<2>::compute_index(const unsigned int n,
                                           unsigned int      (&index)[2]) const;
template <>
void
TensorProductPolynomials<3>::compute_index(const unsigned int n,
                                           unsigned int      (&index)[3]) const;


/* ---------------- template and inline functions ---------- */

template <int dim>
inline
unsigned int
TensorProductPolynomials<dim>::
x_to_the_dim (const unsigned int x)
{
  unsigned int y = 1;
  for (unsigned int d=0; d<dim; ++d)
    y *= x;
  return y;
}



template <int dim>
template <class Pol>
TensorProductPolynomials<dim>::
TensorProductPolynomials(const std::vector<Pol> &pols)
		:
		polynomials (pols.begin(), pols.end()),
		n_tensor_pols(x_to_the_dim(pols.size()))
{}



template <>
void
AnisotropicPolynomials<1>::compute_index(const unsigned int n,
                                         unsigned int      (&index)[1]) const;
template <>
void
AnisotropicPolynomials<2>::compute_index(const unsigned int n,
                                         unsigned int      (&index)[2]) const;
template <>
void
AnisotropicPolynomials<3>::compute_index(const unsigned int n,
                                         unsigned int      (&index)[3]) const;

/// @endif

#endif
