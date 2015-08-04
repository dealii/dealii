// ---------------------------------------------------------------------
// $Id: fe_q.cc 30037 2013-07-18 16:55:40Z maier $
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/base/polynomials_bernstein.h>

#include <vector>
#include <sstream>

DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
FE_Bernstein<dim,spacedim>::FE_Bernstein (const unsigned int degree)
  :
  FE_Q_Base<TensorProductPolynomials<dim>, dim, spacedim> (
   //TensorProductPolynomials<dim>(dealii::generate_complete_bernstein_basis<double>(degree)),
   this->renumber_bases(degree),
   FiniteElementData<dim>(this->get_dpo_vector(degree),
                          1, degree,
                          FiniteElementData<dim>::H1),
   std::vector<bool> (1, false))
{}


template <int dim, int spacedim>
std::string
FE_Bernstein<dim,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_Bernstein<" << dim << ">(" << this->degree << ")";
  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_Bernstein<dim,spacedim>::clone() const
{
  return new FE_Bernstein<dim,spacedim>(*this);
}


/**
 * Only the assertion differs from the same function in FE_Q_Base!!
 */
template <int dim, int spacedim>
std::vector<unsigned int>
FE_Bernstein<dim,spacedim>::get_dpo_vector(const unsigned int deg)
{
  AssertThrow(deg>0,ExcMessage("FE_Bernstein needs to be of degree > 0."));
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}


template <int dim, int spacedim>
TensorProductPolynomials<dim>
FE_Bernstein<dim, spacedim>::renumber_bases(const unsigned int deg)
{
  TensorProductPolynomials<dim> tpp(dealii::generate_complete_bernstein_basis<double>(deg));
  std::vector<unsigned int> renumber(Utilities::fixed_power<dim>(deg+1));
  const FiniteElementData<dim> fe(this->get_dpo_vector(deg),1,
                                  deg);
  FETools::hierarchic_to_lexicographic_numbering (fe, renumber);
  tpp.set_numbering(renumber);
  return tpp;
}


// explicit instantiations
#include "fe_bernstein.inst"
//#include "../base/polynomials_bernstein.inst"

DEAL_II_NAMESPACE_CLOSE
