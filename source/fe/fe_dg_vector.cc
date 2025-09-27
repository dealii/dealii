// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/polynomials_abf.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/polynomials_nedelec.h>
#include <deal.II/base/polynomials_raviart_thomas.h>

#include <deal.II/fe/fe_dg_vector.templates.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
FE_DGNedelec<dim, spacedim>::FE_DGNedelec(const unsigned int p)
  : FE_DGVector<PolynomialsNedelec<dim>, dim, spacedim>(p, {mapping_nedelec})
{}


template <int dim, int spacedim>
std::string
FE_DGNedelec<dim, spacedim>::get_name() const
{
  // note that FETools::get_fe_by_name() depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGNedelec<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree - 1 << ")";

  return namebuf.str();
}


template <int dim, int spacedim>
FE_DGRaviartThomas<dim, spacedim>::FE_DGRaviartThomas(const unsigned int p)
  : FE_DGVector<PolynomialsRaviartThomas<dim>, dim, spacedim>(
      p,
      {mapping_raviart_thomas})
{}


template <int dim, int spacedim>
std::string
FE_DGRaviartThomas<dim, spacedim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGRaviartThomas<" << Utilities::dim_string(dim, spacedim)
          << ">(" << this->degree - 1 << ")";

  return namebuf.str();
}


template <int dim, int spacedim>
FE_DGBDM<dim, spacedim>::FE_DGBDM(const unsigned int p)
  : FE_DGVector<PolynomialsBDM<dim>, dim, spacedim>(p, {mapping_bdm})
{}


template <int dim, int spacedim>
std::string
FE_DGBDM<dim, spacedim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGBDM<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree - 1 << ")";

  return namebuf.str();
}


#include "fe/fe_dg_vector.inst"

DEAL_II_NAMESPACE_CLOSE
