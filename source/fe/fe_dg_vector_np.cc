// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2014 by the deal.II authors
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


#include <deal.II/fe/fe_dg_vector_np.templates.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/polynomials_acfull.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
FE_ACRed<dim, spacedim>::FE_ACRed (const unsigned int p)
  : FE_DGVector_NP<PolynomialsBDM<dim>, dim, spacedim>(p, mapping_bdm)
{}


template <int dim, int spacedim>
std::string
FE_ACRed<dim, spacedim>::get_name () const
{
  // note that the
  // FETools::get_fe_from_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_ACRed<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree-1 << ")";

  return namebuf.str();
}

template <int dim, int spacedim>
FE_ACFull<dim, spacedim>::FE_ACFull (const unsigned int p)
  : FE_DGVector_NP<PolynomialsACFull<dim>, dim, spacedim>(p, mapping_bdm)
{}


template <int dim, int spacedim>
std::string
FE_ACFull<dim, spacedim>::get_name () const
{
  // note that the
  // FETools::get_fe_from_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_ACFull<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree-1 << ")";

  return namebuf.str();
}


#include "fe_dg_vector_np.inst"

DEAL_II_NAMESPACE_CLOSE
