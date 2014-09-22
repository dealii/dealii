// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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


#include <deal.II/base/polynomials_p.h>

DEAL_II_NAMESPACE_OPEN


template <int dim>
PolynomialsP<dim>::PolynomialsP (const unsigned int p)
  :
  PolynomialSpace<dim>(Polynomials::Monomial<double>::generate_complete_basis(p)),
  p(p)
{
  std::vector<unsigned int> index_map(this->n());
  create_polynomial_ordering(index_map);
  this->set_numbering(index_map);
}


template <>
void PolynomialsP<1>::create_polynomial_ordering(
  std::vector<unsigned int> &index_map) const
{
  Assert(index_map.size()==this->n(),
         ExcDimensionMismatch(index_map.size(), this->n()));

  // identity
  for (unsigned int i=0; i<this->n(); ++i)
    index_map[i]=i;
}



const unsigned int imap2[6][21]=
{
  {0},
  {0,1,2},
  {0,1,3,4,2,5},
  {0,1,4,5,2,7,6,8,3,9},
  {0,1,5,6,2,9,7,10,3,12,11,8,13,4,14},
  {0,1,6,7,2,11,8,12,3,15,13,9,16,4,18,14,17,10,19,5,20}
};


template <>
void PolynomialsP<2>::create_polynomial_ordering(
  std::vector<unsigned int> &index_map) const
{
  Assert(index_map.size()==this->n(),
         ExcDimensionMismatch(index_map.size(), this->n()));
  Assert(p<=5, ExcNotImplemented());

  // Given the number i of the
  // polynomial in
  // $1,x,y,xy,x2,y2,...$,
  // index_map[i] gives the number of
  // the polynomial in
  // PolynomialSpace.
  for (unsigned int i=0; i<this->n(); ++i)
    index_map[i]=imap2[p][i];
}


const unsigned int imap3[4][20]=
{
  {0},
  {0,1,2,3},
  {0,1,3,6,4,7,8,2,5,9},
  {0,1,4,10,5,11,13,2,7,16,14,6,12,8,15,17,18,3,9,19}
};

template <>
void PolynomialsP<3>::create_polynomial_ordering(
  std::vector<unsigned int> &index_map) const
{
  Assert(index_map.size()==this->n(),
         ExcDimensionMismatch(index_map.size(), this->n()));
  Assert(p<=3, ExcNotImplemented());

  // Given the number i of the
  // polynomial in
  // $1,x,y,xy,x2,y2,...$,
  // index_map[i] gives the number of
  // the polynomial in
  // PolynomialSpace.
  for (unsigned int i=0; i<this->n(); ++i)
    index_map[i]=imap3[p][i];
}



template class PolynomialsP<1>;
template class PolynomialsP<2>;
template class PolynomialsP<3>;

DEAL_II_NAMESPACE_CLOSE
