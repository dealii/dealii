//-------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------
#ifndef __deal2__polynomials_BDM_h
#define __deal2__polynomials_BDM_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/tensor.h>
#include <base/point.h>
#include <base/polynomial.h>
#include <base/polynomial_space.h>
#include <base/table.h>

#include <vector>


/**
 * @brief The set of BDM polynomials on tensor product cells
 *
 * This class implements the <I>H<SUB>div</SUB></I>-conforming,
 * vector-valued Brezzi-Douglas-Marini polynomials as described in the
 * book by Brezzi and Fortin.
 *
 * Right now, they are implemented in two dimensions only. There, they
 * consist of the complete polynomial space of degree $p$ plus two
 * additional vectors.
 *
 * @author Guido Kanschat, 2003
 */
template <int dim>
class PolynomialsBDM
{
  public:
				     /**
				      * Constructor. Creates all basis
				      * functions for BDM polynomials
				      * of given degree.
				      *
				      * Remark that the degree of a
				      * BDM space is the degree of the
				      * largest complete polynomial
				      * space embedded.
				      *
				      * @arg p: the degree of the
				      * BDM-space
				      */
    PolynomialsBDM (const unsigned int p);

				     /**
				      * Computes the value and the
				      * first and second derivatives
				      * of each BDM
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
    void compute (const Point<dim>            &unit_point,
                  std::vector<Tensor<1,dim> > &values,
                  std::vector<Tensor<2,dim> > &grads,
                  std::vector<Tensor<3,dim> > &grad_grads) const;
    
				     /**
				      * Computes the value of the
				      * <tt>i</tt>th BDM
				      * polynomial at
				      * <tt>unit_point</tt>.
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
    Tensor<1,dim> compute_value (const unsigned int i,
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
    Tensor<2,dim> compute_grad (const unsigned int i,
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
    Tensor<3,dim> compute_grad_grad (const unsigned int i,
                                     const Point<dim> &p) const;

				     /**
				      * Compute the matrix that has as
				      * its entry
				      * <i>a<sub>ij</sub></i> the node
				      * functional <i>i</i> evaluated
				      * for basis function
				      * <i>j</i>. The node functionals
				      * are the standard BDM
				      * interpolation operators.
				      *
				      * The inverse of this matrix can
				      * be used to interpolate node
				      * values to BDM polynomials.
				      */
    void compute_node_matrix (Table<2,double>&) const;
    
				     /**
				      * Returns the number of BDM polynomials.
				      */
    unsigned int n () const;
    
				     /**
				      * Returns the degree of the BDM
				      * space, which is one less than
				      * the highest polynomial degree.
				      */
    unsigned int degree () const;
    
				     /**
				      * Return the number of
				      * polynomials in the space
				      * <TT>BDM(degree)</tt> without
				      * requiring to build an object
				      * of PolynomialsBDM. This is
				      * required by the FiniteElement
				      * classes.
				      */
    static unsigned int compute_n_pols(unsigned int degree);
    
  private:
				     /**
				      * An object representing the
				      * polynomial space used
				      * here. The constructor fills
				      * this with the monomial basis.
				      */
    const PolynomialSpace<dim> polynomial_space;

				     /**
				      * Storage for monomials
				      */
    std::vector<Polynomials::Polynomial<double> > monomials;
    
				     /**
				      * Number of BDM
				      * polynomials.
				      */
    unsigned int n_pols;
    
				     /**
				      * Auxiliary memory.
				      */
    mutable std::vector<double> p_values;
				     /**
				      * Auxiliary memory.
				      */
    mutable std::vector<Tensor<1,dim> > p_grads;
				     /**
				      * Auxiliary memory.
				      */
    mutable std::vector<Tensor<2,dim> > p_grad_grads;
};



template <int dim>
inline unsigned int
PolynomialsBDM<dim>::n() const
{
  return n_pols;
}

template <int dim>
inline unsigned int
PolynomialsBDM<dim>::degree() const
{
  return polynomial_space.degree() - 1;
}

#endif
