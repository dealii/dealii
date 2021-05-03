// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/polynomials_p.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/polynomials_rannacher_turek.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/tensor_product_polynomials_bubbles.h>
#include <deal.II/base/tensor_product_polynomials_const.h>

#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_poly.templates.h>
#include <deal.II/fe/fe_values.h>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
FE_Poly<dim, spacedim>::FE_Poly(const FE_Poly &fe)
  : FiniteElement<dim, spacedim>(fe)
  , poly_space(fe.poly_space->clone())
{}



#include "fe_poly.inst"

DEAL_II_NAMESPACE_CLOSE
