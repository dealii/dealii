// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2021 by the deal.II authors
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

#include <deal.II/base/array_view.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>

#include <memory>
#include <numeric>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
MappingQ<dim, spacedim>::MappingQ(const unsigned int degree)
  : MappingQGeneric<dim, spacedim>(degree)
{}



template <int dim, int spacedim>
MappingQ<dim, spacedim>::MappingQ(const unsigned int degree, const bool)
  : MappingQGeneric<dim, spacedim>(degree)
{}



template <int dim, int spacedim>
MappingQ<dim, spacedim>::MappingQ(const MappingQ<dim, spacedim> &mapping)
  : MappingQGeneric<dim, spacedim>(mapping)
{}


// explicit instantiations
#include "mapping_q.inst"


DEAL_II_NAMESPACE_CLOSE
