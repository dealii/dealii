// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


#include <deal.II/fe/fe_dg_vector.templates.h>
#include <deal.II/base/polynomials_abf.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/polynomials_nedelec.h>
#include <deal.II/base/polynomials_raviart_thomas.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
FE_DGNedelec<dim, spacedim>::FE_DGNedelec (const unsigned int p)
  : FE_DGVector<PolynomialsNedelec<dim>, dim, spacedim>(p, mapping_nedelec)
{}


template <int dim, int spacedim>
std::string
FE_DGNedelec<dim, spacedim>::get_name () const
{
  // note that the
  // FETools::get_fe_from_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGNedelec<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree-1 << ")";

  return namebuf.str();
}


template <int dim, int spacedim>
FE_DGRaviartThomas<dim, spacedim>::FE_DGRaviartThomas (const unsigned int p)
  : FE_DGVector<PolynomialsRaviartThomas<dim>, dim, spacedim>(p, mapping_raviart_thomas)
{}


template <int dim, int spacedim>
std::string
FE_DGRaviartThomas<dim, spacedim>::get_name () const
{
  // note that the
  // FETools::get_fe_from_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGRaviartThomas<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree-1 << ")";

  return namebuf.str();
}


template <int dim, int spacedim>
FE_DGBDM<dim, spacedim>::FE_DGBDM (const unsigned int p)
  : FE_DGVector<PolynomialsBDM<dim>, dim, spacedim>(p, mapping_bdm)
{}


template <int dim, int spacedim>
std::string
FE_DGBDM<dim, spacedim>::get_name () const
{
  // note that the
  // FETools::get_fe_from_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGBDM<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree-1 << ")";

  return namebuf.str();
}


#include "fe_dg_vector.inst"

DEAL_II_NAMESPACE_CLOSE

