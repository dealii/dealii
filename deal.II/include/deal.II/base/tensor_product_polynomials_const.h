//---------------------------------------------------------------------------
//    $Id: tensor_product_polynomials.h 27628 2012-11-20 22:49:26Z heister $
//
//    Copyright (C) 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__tensor_product_polynomials_const_h
#define __deal2__tensor_product_polynomials_const_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Polynomials
 * @{
 */

/**
 * Tensor product of given polynomials.
 *
 * Given a vector of <i>n</i> one-dimensional polynomials
 * <i>P<sub>1</sub></i> to <i>P<sub>n</sub></i>, this class generates
 * <i>n<sup>dim</sup></i> polynomials of the form
 * <i>Q<sub>ijk</sub>(x,y,z) =
 * P<sub>i</sub>(x)P<sub>j</sub>(y)P<sub>k</sub>(z)</i>and a locally
 * constant function. If the base polynomials are mutually orthogonal
 * on the interval [-1,1] or [0,1], then the tensor product
 * polynomials are orthogonal on [-1,1]<sup>dim</sup> or [0,1]
 * <sup>dim</sup>, respectively.
 *
 * Indexing is as follows: the order of dim-dimensional polynomials is
 * x-coordinates running fastest, then y-coordinate, etc. The first
 * few polynomials are thus <i>P<sub>1</sub>(x)P<sub>1</sub>(y),
 * P<sub>2</sub>(x)P<sub>1</sub>(y), P<sub>3</sub>(x)P<sub>1</sub>(y),
 * ..., P<sub>1</sub>(x)P<sub>2</sub>(y),
 * P<sub>2</sub>(x)P<sub>2</sub>(y), P<sub>3</sub>(x)P<sub>2</sub>(y),
 * ...</i> and likewise in 3d. The locally constant function has the last index
 *
 * The output_indices() function prints the ordering of the
 * dim-dimensional polynomials, i.e. for each polynomial in the
 * polynomial space it gives the indices i,j,k of the one-dimensional
 * polynomials in x,y and z direction. The ordering of the
 * dim-dimensional polynomials can be changed by using the
 * set_numbering() function.
 *
 * @author Ralf Hartmann, 2000, 2004, Guido Kanschat, 2000, Wolfgang Bangerth 2003
 */
template <int dim>
class TensorProductPolynomialsConst
{
public:
  /**
   * Access to the dimension of
   * this object, for checking and
   * automatic setting of dimension
   * in other classes.
   */
  static const unsigned int dimension = dim;

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
  TensorProductPolynomialsConst (const std::vector<Pol> &pols);

  /**
   * Prints the list of the indices
   * to <tt>out</tt>.
   */
  void output_indices(std::ostream &out) const;

  /**
   * Sets the ordering of the
   * polynomials. Requires
   * <tt>renumber.size()==n()</tt>.
   * Stores a copy of
   * <tt>renumber</tt>.
   */
  void set_numbering(const std::vector<unsigned int> &renumber);

  /**
   * Gives read access to the
   * renumber vector.
   */
  const std::vector<unsigned int> &get_numbering() const;

  /**
   * Gives read access to the
   * inverse renumber vector.
   */
  const std::vector<unsigned int> &get_numbering_inverse() const;

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
   * product polynomials plus the constant
   * function. For <i>n</i> 1d polynomials
   * this is <i>n<sup>dim</sup>+1</i>.
   */
  unsigned int n () const;


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
   * Index map for reordering the
   * polynomials.
   */
  std::vector<unsigned int> index_map;

  /**
   * Index map for reordering the
   * polynomials.
   */
  std::vector<unsigned int> index_map_inverse;

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
  // fix to avoid compiler warnings about zero
  // length arrays
  void compute_index (const unsigned int i,
                      unsigned int       (&indices)[(dim>0?dim:1)]) const;

  /**
   * Computes
   * <i>x<sup>dim</sup></i> for
   * unsigned int <i>x</i>. Used in
   * the constructor.
   */
  static
  unsigned int x_to_the_dim (const unsigned int x);
};

#ifndef DOXYGEN

template <int dim>
inline
const std::vector<unsigned int> &
TensorProductPolynomialsConst<dim>::get_numbering() const
{
  return index_map;
}


template <int dim>
inline
const std::vector<unsigned int> &
TensorProductPolynomialsConst<dim>::get_numbering_inverse() const
{
  return index_map_inverse;
}


#endif // DOXYGEN

#ifndef DOXYGEN

/* -------------- declaration of explicit specializations --- */

template <>
void
TensorProductPolynomialsConst<1>::compute_index(const unsigned int n,
                                                unsigned int      (&index)[1]) const;
template <>
void
TensorProductPolynomialsConst<2>::compute_index(const unsigned int n,
                                                unsigned int      (&index)[2]) const;
template <>
void
TensorProductPolynomialsConst<3>::compute_index(const unsigned int n,
                                                unsigned int      (&index)[3]) const;


/* ---------------- template and inline functions ---------- */

template <int dim>
inline
unsigned int
TensorProductPolynomialsConst<dim>::
x_to_the_dim (const unsigned int x)
{
  unsigned int y = 1;
  for (int d=0; d<dim; ++d)
    y *= x;
  return y;
}



template <int dim>
template <class Pol>
TensorProductPolynomialsConst<dim>::
TensorProductPolynomialsConst(const std::vector<Pol> &pols)
  :
  polynomials (pols.begin(), pols.end()),
  n_tensor_pols(x_to_the_dim(pols.size())+1),
  index_map(n_tensor_pols),
  index_map_inverse(n_tensor_pols)
{
  // per default set this index map
  // to identity. This map can be
  // changed by the user through the
  // set_numbering() function
  for (unsigned int i=0; i<n_tensor_pols; ++i)
    {
      index_map[i]=i;
      index_map_inverse[i]=i;
    }
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
