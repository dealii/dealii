// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__abf_scalar_polynomials_h
#define __deal2__abf_scalar_polynomials_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/utilities.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Polynomials
 * @{
 */

/**
 *
 * @author Zhen Tao, 2015
 */
template <int dim, typename POLY=Polynomials::Polynomial<double> >
class ABFScalarPolynomials
{
public:
  /**
   * Access to the dimension of this object, for checking and automatic
   * setting of dimension in other classes.
   */
  static const unsigned int dimension = dim;

  /**
   * Constructor. <tt>pols</tt> is a vector of objects that should be derived
   * or otherwise convertible to one-dimensional polynomial objects of type @p
   * POLY (template argument of class). It will be copied element by element
   * into a private variable.

   * Note: Polynomials have to be Legendre, not LagrangeEquidistant
   */
  template <class Pol>
  ABFScalarPolynomials (const std::vector<Pol> &pols);

  /**
   * Prints the list of the indices to <tt>out</tt>.
   */
  void output_indices(std::ostream &out) const;

  /**
   * Sets the ordering of the polynomials. Requires
   * <tt>renumber.size()==n()</tt>.  Stores a copy of <tt>renumber</tt>.
   */
  void set_numbering(const std::vector<unsigned int> &renumber);

  /**
   * Gives read access to the renumber vector.
   */
  const std::vector<unsigned int> &get_numbering() const;

  /**
   * Gives read access to the inverse renumber vector.
   */
  const std::vector<unsigned int> &get_numbering_inverse() const;

  /**
   * Computes the value and the first and second derivatives of each tensor
   * product polynomial at <tt>unit_point</tt>.
   *
   * The size of the vectors must either be equal 0 or equal n(). In the first
   * case, the function will not compute these values.
   *
   * If you need values or derivatives of all tensor product polynomials then
   * use this function, rather than using any of the compute_value(),
   * compute_grad() or compute_grad_grad() functions, see below, in a loop
   * over all tensor product polynomials.
   */
  void compute (const Point<dim>            &unit_point,
                std::vector<double>         &values,
                std::vector<Tensor<1,dim> > &grads,
                std::vector<Tensor<2,dim> > &grad_grads) const;

  /**
   * Computes the value of the <tt>i</tt>th tensor product polynomial at
   * <tt>unit_point</tt>. Here <tt>i</tt> is given in tensor product
   * numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each point value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the compute() function with
   * <tt>values.size()==</tt>n() to get the point values of all tensor
   * polynomials all at once and in a much more efficient way.
   */
  double compute_value (const unsigned int i,
                        const Point<dim> &p) const;

  /**
   * Computes the grad of the <tt>i</tt>th tensor product polynomial at
   * <tt>unit_point</tt>. Here <tt>i</tt> is given in tensor product
   * numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the compute() function, see above, with
   * <tt>grads.size()==</tt>n() to get the point value of all tensor
   * polynomials all at once and in a much more efficient way.
   */
  Tensor<1,dim> compute_grad (const unsigned int i,
                              const Point<dim> &p) const;

  /**
   * Computes the second derivative (grad_grad) of the <tt>i</tt>th tensor
   * product polynomial at <tt>unit_point</tt>. Here <tt>i</tt> is given in
   * tensor product numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the compute() function, see above, with
   * <tt>grad_grads.size()==</tt>n() to get the point value of all tensor
   * polynomials all at once and in a much more efficient way.
   */
  Tensor<2,dim> compute_grad_grad (const unsigned int i,
                                   const Point<dim> &p) const;

  /**
   * Returns the number of tensor product polynomials. For <i>n</i> 1d
   * polynomials this is <i>n<sup>dim</sup></i>.
   */
  unsigned int n () const;


protected:
  /**
   * Copy of the vector <tt>pols</tt> of polynomials given to the constructor.
   */
  std::vector<POLY> polynomials;

  /**
   * Number of tensor product polynomials. See n().
   */
  unsigned int n_tensor_pols;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map;

  /**
   * Index map for reordering the polynomials.
   */
  std::vector<unsigned int> index_map_inverse;

  /**
   * Each tensor product polynomial <i>i</i> is a product of one-dimensional
   * polynomials in each space direction. Compute the indices of these one-
   * dimensional polynomials for each space direction, given the index
   * <i>i</i>.
   */
  // fix to avoid compiler warnings about zero length arrays
  void compute_index (const unsigned int i,
                      unsigned int       (&indices)[(dim>0?dim:1)]) const;

  static
  unsigned int
  get_n_tensor_pols (const std::vector<POLY> &pols);

};

#ifndef DOXYGEN


/* ---------------- template and inline functions ---------- */


template <int dim, typename POLY>
template <class Pol>
inline
ABFScalarPolynomials<dim,POLY>::
ABFScalarPolynomials(const std::vector<Pol> &pols)
  :
  polynomials (pols.begin(), pols.end()),
  n_tensor_pols(get_n_tensor_pols(pols)),
  index_map(n_tensor_pols),
  index_map_inverse(n_tensor_pols)
{
  // per default set this index map to identity. This map can be changed by
  // the user through the set_numbering() function
  for (unsigned int i=0; i<n_tensor_pols; ++i)
    {
      index_map[i]=i;
      index_map_inverse[i]=i;
    }
}



template <int dim, typename POLY>
inline
unsigned int
ABFScalarPolynomials<dim,POLY>::n() const
{
  if (dim == 0)
    return numbers::invalid_unsigned_int;
  else
    return n_tensor_pols;
}



template <int dim, typename POLY>
inline
const std::vector<unsigned int> &
ABFScalarPolynomials<dim,POLY>::get_numbering() const
{
  return index_map;
}


template <int dim, typename POLY>
inline
const std::vector<unsigned int> &
ABFScalarPolynomials<dim,POLY>::get_numbering_inverse() const
{
  return index_map_inverse;
}



#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
