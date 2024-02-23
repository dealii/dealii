// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

#ifndef DOXYGEN

template <int dim, int spacedim>
FE_Poly<dim, spacedim>::FE_Poly(const FE_Poly &fe)
  : FiniteElement<dim, spacedim>(fe)
  , poly_space(fe.poly_space->clone())
{}

#endif

#include "fe_poly.inst"

DEAL_II_NAMESPACE_CLOSE
