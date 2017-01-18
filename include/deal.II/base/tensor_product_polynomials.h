// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2016 by the deal.II authors
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

#ifndef dealii__tensor_product_polynomials_h
#define dealii__tensor_product_polynomials_h


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
 * Tensor product of given polynomials.
 *
 * Given a vector of <i>n</i> one-dimensional polynomials <i>P<sub>1</sub></i>
 * to <i>P<sub>n</sub></i>, this class generates <i>n<sup>dim</sup></i>
 * polynomials of the form <i>Q<sub>ijk</sub>(x,y,z) =
 * P<sub>i</sub>(x)P<sub>j</sub>(y)P<sub>k</sub>(z)</i>. If the base
 * polynomials are mutually orthogonal on the interval [-1,1] or [0,1], then
 * the tensor product polynomials are orthogonal on [-1,1]<sup>dim</sup> or
 * [0,1]<sup>dim</sup>, respectively.
 *
 * Indexing is as follows: the order of dim-dimensional polynomials is
 * x-coordinates running fastest, then y-coordinate, etc. The first few
 * polynomials are thus <i>P<sub>1</sub>(x)P<sub>1</sub>(y),
 * P<sub>2</sub>(x)P<sub>1</sub>(y), P<sub>3</sub>(x)P<sub>1</sub>(y), ...,
 * P<sub>1</sub>(x)P<sub>2</sub>(y), P<sub>2</sub>(x)P<sub>2</sub>(y),
 * P<sub>3</sub>(x)P<sub>2</sub>(y), ...</i> and likewise in 3d.
 *
 * The output_indices() function prints the ordering of the dim-dimensional
 * polynomials, i.e. for each polynomial in the polynomial space it gives the
 * indices i,j,k of the one-dimensional polynomials in x,y and z direction.
 * The ordering of the dim-dimensional polynomials can be changed by using the
 * set_numbering() function.
 *
 * @author Ralf Hartmann, 2000, 2004, Guido Kanschat, 2000, Wolfgang Bangerth
 * 2003
 */
template <int dim, typename PolynomialType=Polynomials::Polynomial<double> >
class TensorProductPolynomials
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
   * PolynomialType (template argument of class). It will be copied element by
   * element into a private variable.
   */
  template <class Pol>
  TensorProductPolynomials (const std::vector<Pol> &pols);

  /**
   * Prints the list of the indices to <tt>out</tt>.
   */
  void output_indices(std::ostream &out) const;

  /**
   * Set the ordering of the polynomials. Requires
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
   * Compute the value and the first and second derivatives of each tensor
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
                std::vector<Tensor<2,dim> > &grad_grads,
                std::vector<Tensor<3,dim> > &third_derivatives,
                std::vector<Tensor<4,dim> > &fourth_derivatives) const;

  /**
   * Compute the value of the <tt>i</tt>th tensor product polynomial at
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
   * Compute the <tt>order</tt>th derivative of the <tt>i</tt>th tensor
   * product polynomial at <tt>unit_point</tt>. Here <tt>i</tt> is given in
   * tensor product numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the compute() function, see above, with the
   * size of the appropriate parameter set to n() to get the point value of
   * all tensor polynomials all at once and in a much more efficient way.
   *
   * @tparam order The derivative order.
   */
  template <int order>
  Tensor<order,dim> compute_derivative (const unsigned int i,
                                        const Point<dim> &p) const;

  /**
   * Compute the grad of the <tt>i</tt>th tensor product polynomial at
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
   * Compute the second derivative (grad_grad) of the <tt>i</tt>th tensor
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
   * Return the number of tensor product polynomials. For <i>n</i> 1d
   * polynomials this is <i>n<sup>dim</sup></i>.
   */
  unsigned int n () const;


protected:
  /**
   * Copy of the vector <tt>pols</tt> of polynomials given to the constructor.
   */
  std::vector<PolynomialType> polynomials;

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
};



/**
 * Anisotropic tensor product of given polynomials.
 *
 * Given one-dimensional polynomials <tt>Px1</tt>, <tt>Px2</tt>, ... in
 * x-direction, <tt>Py1</tt>, <tt>Py2</tt>, ... in y-direction, and so on,
 * this class generates polynomials of the form  <i>Q<sub>ijk</sub>(x,y,z) =
 * Pxi(x)Pyj(y)Pzk(z)</i>. If the base polynomials are mutually orthogonal on
 * the interval $[-1,1]$ or $[0,d]$, then the tensor product polynomials are
 * orthogonal on $[-1,1]^d$ or $[0,1]^d$, respectively.
 *
 * Indexing is as follows: the order of dim-dimensional polynomials is
 * x-coordinates running fastest, then y-coordinate, etc. The first few
 * polynomials are thus <tt>Px1(x)Py1(y)</tt>, <tt>Px2(x)Py1(y)</tt>,
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
   * Constructor. <tt>pols</tt> is a table of one-dimensional polynomials. The
   * number of rows in this table should be equal to the space dimension, with
   * the elements of each row giving the polynomials that shall be used in
   * this particular coordinate direction. These polynomials may vary between
   * coordinates, as well as their number.
   */
  AnisotropicPolynomials (const std::vector<std::vector<Polynomials::Polynomial<double> > > &pols);

  /**
   * Compute the value and the first and second derivatives of each tensor
   * product polynomial at <tt>unit_point</tt>.
   *
   * The size of the vectors must either be equal <tt>0</tt> or equal
   * <tt>n_tensor_pols</tt>.  In the first case, the function will not compute
   * these values.
   *
   * If you need values or derivatives of all tensor product polynomials then
   * use this function, rather than using any of the <tt>compute_value</tt>,
   * <tt>compute_grad</tt> or <tt>compute_grad_grad</tt> functions, see below,
   * in a loop over all tensor product polynomials.
   */
  void compute (const Point<dim>            &unit_point,
                std::vector<double>         &values,
                std::vector<Tensor<1,dim> > &grads,
                std::vector<Tensor<2,dim> > &grad_grads,
                std::vector<Tensor<3,dim> > &third_derivatives,
                std::vector<Tensor<4,dim> > &fourth_derivatives) const;

  /**
   * Compute the value of the <tt>i</tt>th tensor product polynomial at
   * <tt>unit_point</tt>. Here <tt>i</tt> is given in tensor product
   * numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each point value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the <tt>compute</tt> function, see above,
   * with <tt>values.size()==n_tensor_pols</tt> to get the point values of all
   * tensor polynomials all at once and in a much more efficient way.
   */
  double compute_value (const unsigned int i,
                        const Point<dim> &p) const;

  /**
   * Compute the <tt>order</tt>th derivative of the <tt>i</tt>th tensor
   * product polynomial at <tt>unit_point</tt>. Here <tt>i</tt> is given in
   * tensor product numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the compute() function, see above, with the
   * size of the appropriate parameter set to n() to get the point value of
   * all tensor polynomials all at once and in a much more efficient way.
   *
   * @tparam order The derivative order.
   */
  template <int order>
  Tensor<order,dim> compute_derivative (const unsigned int i,
                                        const Point<dim> &p) const;

  /**
   * Compute the grad of the <tt>i</tt>th tensor product polynomial at
   * <tt>unit_point</tt>. Here <tt>i</tt> is given in tensor product
   * numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the <tt>compute</tt> function, see above,
   * with <tt>grads.size()==n_tensor_pols</tt> to get the point value of all
   * tensor polynomials all at once and in a much more efficient way.
   */
  Tensor<1,dim> compute_grad (const unsigned int i,
                              const Point<dim> &p) const;

  /**
   * Compute the second derivative (grad_grad) of the <tt>i</tt>th tensor
   * product polynomial at <tt>unit_point</tt>. Here <tt>i</tt> is given in
   * tensor product numbering.
   *
   * Note, that using this function within a loop over all tensor product
   * polynomials is not efficient, because then each derivative value of the
   * underlying (one-dimensional) polynomials is (unnecessarily) computed
   * several times.  Instead use the <tt>compute</tt> function, see above,
   * with <tt>grad_grads.size()==n_tensor_pols</tt> to get the point value of
   * all tensor polynomials all at once and in a much more efficient way.
   */
  Tensor<2,dim> compute_grad_grad (const unsigned int i,
                                   const Point<dim> &p) const;

  /**
   * Return the number of tensor product polynomials. It is the product of
   * the number of polynomials in each coordinate direction.
   */
  unsigned int n () const;

private:
  /**
   * Copy of the vector <tt>pols</tt> of polynomials given to the constructor.
   */
  const std::vector<std::vector<Polynomials::Polynomial<double> > > polynomials;

  /**
   * Number of tensor product polynomials. This is <tt>Nx*Ny*Nz</tt>, or with
   * terms dropped if the number of space dimensions is less than 3.
   */
  const unsigned int n_tensor_pols;

  /**
   * Each tensor product polynomial @þ{i} is a product of one-dimensional
   * polynomials in each space direction. Compute the indices of these one-
   * dimensional polynomials for each space direction, given the index
   * <tt>i</tt>.
   */
  void compute_index (const unsigned int i,
                      unsigned int       (&indices)[dim]) const;

  /**
   * Given the input to the constructor, compute <tt>n_tensor_pols</tt>.
   */
  static
  unsigned int
  get_n_tensor_pols (const std::vector<std::vector<Polynomials::Polynomial<double> > > &pols);
};

/** @} */

#ifndef DOXYGEN


/* ---------------- template and inline functions ---------- */


template <int dim, typename PolynomialType>
template <class Pol>
inline
TensorProductPolynomials<dim,PolynomialType>::
TensorProductPolynomials(const std::vector<Pol> &pols)
  :
  polynomials (pols.begin(), pols.end()),
  n_tensor_pols(Utilities::fixed_power<dim>(pols.size())),
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



template <int dim, typename PolynomialType>
inline
unsigned int
TensorProductPolynomials<dim,PolynomialType>::n() const
{
  if (dim == 0)
    return numbers::invalid_unsigned_int;
  else
    return n_tensor_pols;
}



template <int dim, typename PolynomialType>
inline
const std::vector<unsigned int> &
TensorProductPolynomials<dim,PolynomialType>::get_numbering() const
{
  return index_map;
}


template <int dim, typename PolynomialType>
inline
const std::vector<unsigned int> &
TensorProductPolynomials<dim,PolynomialType>::get_numbering_inverse() const
{
  return index_map_inverse;
}

template <int dim, typename PolynomialType>
template <int order>
Tensor<order,dim>
TensorProductPolynomials<dim,PolynomialType>::compute_derivative
(const unsigned int  i,
 const Point<dim>   &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);

  double v [dim][5];
  {
    std::vector<double> tmp (5);
    for (unsigned int d=0; d<dim; ++d)
      {
        polynomials[indices[d]].value (p(d), tmp);
        v[d][0] = tmp[0];
        v[d][1] = tmp[1];
        v[d][2] = tmp[2];
        v[d][3] = tmp[3];
        v[d][4] = tmp[4];
      }
  }

  Tensor<order,dim> derivative;
  switch (order)
    {
    case 1:
    {
      Tensor<1,dim> &derivative_1 = *reinterpret_cast<Tensor<1,dim>*>(&derivative);
      for (unsigned int d=0; d<dim; ++d)
        {
          derivative_1[d] = 1.;
          for (unsigned int x=0; x<dim; ++x)
            {
              unsigned int x_order=0;
              if (d==x) ++x_order;

              derivative_1[d] *= v[x][x_order];
            }
        }

      return derivative;
    }
    case 2:
    {
      Tensor<2,dim> &derivative_2 = *reinterpret_cast<Tensor<2,dim>*>(&derivative);
      for (unsigned int d1=0; d1<dim; ++d1)
        for (unsigned int d2=0; d2<dim; ++d2)
          {
            derivative_2[d1][d2] = 1.;
            for (unsigned int x=0; x<dim; ++x)
              {
                unsigned int x_order=0;
                if (d1==x) ++x_order;
                if (d2==x) ++x_order;

                derivative_2[d1][d2] *= v[x][x_order];
              }
          }

      return derivative;
    }
    case 3:
    {
      Tensor<3,dim> &derivative_3 = *reinterpret_cast<Tensor<3,dim>*>(&derivative);
      for (unsigned int d1=0; d1<dim; ++d1)
        for (unsigned int d2=0; d2<dim; ++d2)
          for (unsigned int d3=0; d3<dim; ++d3)
            {
              derivative_3[d1][d2][d3] = 1.;
              for (unsigned int x=0; x<dim; ++x)
                {
                  unsigned int x_order=0;
                  if (d1==x) ++x_order;
                  if (d2==x) ++x_order;
                  if (d3==x) ++x_order;

                  derivative_3[d1][d2][d3] *= v[x][x_order];
                }
            }

      return derivative;
    }
    case 4:
    {
      Tensor<4,dim> &derivative_4 = *reinterpret_cast<Tensor<4,dim>*>(&derivative);
      for (unsigned int d1=0; d1<dim; ++d1)
        for (unsigned int d2=0; d2<dim; ++d2)
          for (unsigned int d3=0; d3<dim; ++d3)
            for (unsigned int d4=0; d4<dim; ++d4)
              {
                derivative_4[d1][d2][d3][d4] = 1.;
                for (unsigned int x=0; x<dim; ++x)
                  {
                    unsigned int x_order=0;
                    if (d1==x) ++x_order;
                    if (d2==x) ++x_order;
                    if (d3==x) ++x_order;
                    if (d4==x) ++x_order;

                    derivative_4[d1][d2][d3][d4] *= v[x][x_order];
                  }
              }

      return derivative;
    }
    default:
    {
      Assert (false, ExcNotImplemented());
      return derivative;
    }
    }
}

template <int dim>
template <int order>
Tensor<order,dim>
AnisotropicPolynomials<dim>::compute_derivative (const unsigned int i,
                                                 const Point<dim> &p) const
{
  unsigned int indices[dim];
  compute_index (i, indices);

  std::vector<std::vector<double> > v(dim, std::vector<double> (order+1));
  for (unsigned int d=0; d<dim; ++d)
    polynomials[d][indices[d]].value(p(d), v[d]);

  Tensor<order,dim> derivative;
  switch (order)
    {
    case 1:
    {
      Tensor<1,dim> &derivative_1 = *reinterpret_cast<Tensor<1,dim>*>(&derivative);
      for (unsigned int d=0; d<dim; ++d)
        {
          derivative_1[d] = 1.;
          for (unsigned int x=0; x<dim; ++x)
            {
              unsigned int x_order=0;
              if (d==x) ++x_order;

              derivative_1[d] *= v[x][x_order];
            }
        }

      return derivative;
    }
    case 2:
    {
      Tensor<2,dim> &derivative_2 = *reinterpret_cast<Tensor<2,dim>*>(&derivative);
      for (unsigned int d1=0; d1<dim; ++d1)
        for (unsigned int d2=0; d2<dim; ++d2)
          {
            derivative_2[d1][d2] = 1.;
            for (unsigned int x=0; x<dim; ++x)
              {
                unsigned int x_order=0;
                if (d1==x) ++x_order;
                if (d2==x) ++x_order;

                derivative_2[d1][d2] *= v[x][x_order];
              }
          }

      return derivative;
    }
    case 3:
    {
      Tensor<3,dim> &derivative_3 = *reinterpret_cast<Tensor<3,dim>*>(&derivative);
      for (unsigned int d1=0; d1<dim; ++d1)
        for (unsigned int d2=0; d2<dim; ++d2)
          for (unsigned int d3=0; d3<dim; ++d3)
            {
              derivative_3[d1][d2][d3] = 1.;
              for (unsigned int x=0; x<dim; ++x)
                {
                  unsigned int x_order=0;
                  if (d1==x) ++x_order;
                  if (d2==x) ++x_order;
                  if (d3==x) ++x_order;

                  derivative_3[d1][d2][d3] *= v[x][x_order];
                }
            }

      return derivative;
    }
    case 4:
    {
      Tensor<4,dim> &derivative_4 = *reinterpret_cast<Tensor<4,dim>*>(&derivative);
      for (unsigned int d1=0; d1<dim; ++d1)
        for (unsigned int d2=0; d2<dim; ++d2)
          for (unsigned int d3=0; d3<dim; ++d3)
            for (unsigned int d4=0; d4<dim; ++d4)
              {
                derivative_4[d1][d2][d3][d4] = 1.;
                for (unsigned int x=0; x<dim; ++x)
                  {
                    unsigned int x_order=0;
                    if (d1==x) ++x_order;
                    if (d2==x) ++x_order;
                    if (d3==x) ++x_order;
                    if (d4==x) ++x_order;

                    derivative_4[d1][d2][d3][d4] *= v[x][x_order];
                  }
              }

      return derivative;
    }
    default:
    {
      Assert (false, ExcNotImplemented());
      return derivative;
    }
    }
}



#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
