//----------------------  tensor_product_polynomials.h  -------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
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
 * Given a vector of @{n} one-dimensional polynomials @{P1} to @{Pn},
 * this class generates @p{n} to the power of @p{dim} polynomials of
 * the form @p{ Qijk(x,y,z) = Pi(x)Pj(y)Pk(z)}. If the base
 * polynomials are mutually orthogonal on the interval $[-1,1]$ or
 * $[0,d], then the tensor product polynomials are orthogonal on
 * $[-1,1]^d$ or $[0,1]^d$, respectively.
 *
 * Indexing is as following: the order of dim-dimensional polynomials
 * is x-coordinates running fastest, then y-coordinate, etc. The first
 * few polynomials are thus @p{P1(x)P1(y)}, @p{P2(x)P1(y)},
 * @p{P3(x)P1(y)}, ..., @p{P1(x)P2(y)}, @p{P2(x)P2(y)},
 * @p{P3(x)P2(y)}, ..., and likewise in 3d.
 * 
 * @author Ralf Hartmann, Guido Kanschat, 2000, Wolfgang Bangerth 2003
 */
template <int dim>
class TensorProductPolynomials
{
  public:
				     /**
				      * Constructor. @p{pols} is a
				      * vector of objects that should
				      * be derived or otherwise
				      * convertible to one-dimensional
				      * polynomial objects and will be
				      * copied into the member
				      * variable @p{polynomials}.
				      */
    template <class Pol>
    TensorProductPolynomials (const std::vector<Pol> &pols);

				     /**
				      * Computes the value and the
				      * first and second derivatives
				      * of each tensor product
				      * polynomial at @p{unit_point}.
				      *
				      * The size of the vectors must
				      * either be equal @p{0} or equal
				      * @p{n_tensor_pols}.  In the
				      * first case, the function will
				      * not compute these values.
				      *
				      * If you need values or
				      * derivatives of all tensor
				      * product polynomials then use
				      * this function, rather than
				      * using any of the
				      * @p{compute_value},
				      * @p{compute_grad} or
				      * @p{compute_grad_grad}
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
				      * @p{i}th tensor product
				      * polynomial at
				      * @p{unit_point}. Here @p{i} is
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
				      * the @p{compute} function, see
				      * above, with
				      * @p{values.size()==n_tensor_pols}
				      * to get the point values of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    double compute_value (const unsigned int i,
			  const Point<dim> &p) const;

				     /**
				      * Computes the grad of the
				      * @p{i}th tensor product
				      * polynomial at
				      * @p{unit_point}. Here @p{i} is
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
				      * Instead use the @p{compute}
				      * function, see above, with
				      * @p{grads.size()==n_tensor_pols}
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
				      * @p{i}th tensor product
				      * polynomial at
				      * @p{unit_point}. Here @p{i} is
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
				      * Instead use the @p{compute}
				      * function, see above, with
				      * @p{grad_grads.size()==n_tensor_pols}
				      * to get the point value of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    Tensor<2,dim> compute_grad_grad (const unsigned int i,
                                     const Point<dim> &p) const;

				     /**
				      * Returns the number of tensor
				      * product polynomials. For $n$
				      * 1d polynomials this is $n^dim$.
				      */
    unsigned int n() const;

				     /**
				      * Exception.
				      */
    DeclException3 (ExcDimensionMismatch2,
		    int, int, int,
		    << "Dimension " << arg1 << " not equal to " << arg2 << " nor to " << arg3);

	    
  private:
				     /**
				      * Copy of the vector @p{pols} of
				      * polynomials given to the
				      * constructor.
				      */
    const std::vector<Polynomials::Polynomial<double> > polynomials;

				     /**
				      * Number of tensor product
				      * polynomials. For $n$ 1d
				      * polynomials this is $n^dim$.
				      */
    const unsigned int n_tensor_pols;

                                     /**
                                      * Each tensor product polynomial
                                      * @þ{i} is a product of
                                      * one-dimensional polynomials in
                                      * each space direction. Compute
                                      * the indices of these
                                      * one-dimensional polynomials
                                      * for each space direction,
                                      * given the index @p{i}.
                                      */
    void compute_index (const unsigned int i,
                        unsigned int       (&indices)[dim]) const;
    
				     /**
				      * Computes @p{x} to the power of
				      * @p{dim} for unsigned int @p{x}.
				      * Used in the constructor.
				      */
    static
    unsigned int x_to_the_dim (const unsigned int x);
};




/**
 * Anisotropic tensor product of given polynomials.
 *
 * Given one-dimensional polynomials @{Px1}, @{Px2}, ... in
 * x-direction, @{Py1}, @{Py2}, ... in y-direction, and so on, this
 * class generates polynomials of the form @p{ Qijk(x,y,z) =
 * Pxi(x)Pyj(y)Pzk(z)}. If the base polynomials are mutually
 * orthogonal on the interval $[-1,1]$ or $[0,d], then the tensor
 * product polynomials are orthogonal on $[-1,1]^d$ or $[0,1]^d$,
 * respectively.
 *
 * Indexing is as following: the order of dim-dimensional polynomials
 * is x-coordinates running fastest, then y-coordinate, etc. The first
 * few polynomials are thus @p{Px1(x)Py1(y)}, @p{Px2(x)Py1(y)},
 * @p{Px3(x)Py1(y)}, ..., @p{Px1(x)Py2(y)}, @p{Px2(x)Py2(y)},
 * @p{Px3(x)Py2(y)}, ..., and likewise in 3d.
 * 
 * @author Wolfgang Bangerth 2003
 */
template <int dim>
class AnisotropicPolynomials
{
  public:
				     /**
				      * Constructor. @p{pols} is a
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
				      * polynomial at @p{unit_point}.
				      *
				      * The size of the vectors must
				      * either be equal @p{0} or equal
				      * @p{n_tensor_pols}.  In the
				      * first case, the function will
				      * not compute these values.
				      *
				      * If you need values or
				      * derivatives of all tensor
				      * product polynomials then use
				      * this function, rather than
				      * using any of the
				      * @p{compute_value},
				      * @p{compute_grad} or
				      * @p{compute_grad_grad}
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
				      * @p{i}th tensor product
				      * polynomial at
				      * @p{unit_point}. Here @p{i} is
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
				      * the @p{compute} function, see
				      * above, with
				      * @p{values.size()==n_tensor_pols}
				      * to get the point values of all
				      * tensor polynomials all at once
				      * and in a much more efficient
				      * way.
				      */
    double compute_value (const unsigned int i,
			  const Point<dim> &p) const;

				     /**
				      * Computes the grad of the
				      * @p{i}th tensor product
				      * polynomial at
				      * @p{unit_point}. Here @p{i} is
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
				      * Instead use the @p{compute}
				      * function, see above, with
				      * @p{grads.size()==n_tensor_pols}
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
				      * @p{i}th tensor product
				      * polynomial at
				      * @p{unit_point}. Here @p{i} is
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
				      * Instead use the @p{compute}
				      * function, see above, with
				      * @p{grad_grads.size()==n_tensor_pols}
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
    unsigned int n() const;

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
				      * Copy of the vector @p{pols} of
				      * polynomials given to the
				      * constructor.
				      */
    const std::vector<std::vector<Polynomials::Polynomial<double> > > polynomials;

				     /**
				      * Number of tensor product
				      * polynomials. This is
				      * @p{Nx*Ny*Nz}, or with terms
				      * dropped if the number of space
				      * dimensions is less than 3.
				      */
    const unsigned int n_tensor_pols;

                                     /**
                                      * Each tensor product polynomial
                                      * @þ{i} is a product of
                                      * one-dimensional polynomials in
                                      * each space direction. Compute
                                      * the indices of these
                                      * one-dimensional polynomials
                                      * for each space direction,
                                      * given the index @p{i}.
                                      */
    void compute_index (const unsigned int i,
                        unsigned int       (&indices)[dim]) const;
    
				     /**
				      * Given the input to the
				      * constructor, compute
				      * @p{n_tensor_pols}.
				      */
    static
    unsigned int
    get_n_tensor_pols (const std::vector<std::vector<Polynomials::Polynomial<double> > > &pols);
};




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



#endif
