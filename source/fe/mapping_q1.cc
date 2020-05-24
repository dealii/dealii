// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
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


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>


DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
MappingQ1<dim, spacedim>::MappingQ1()
  : MappingQGeneric<dim, spacedim>(1)
{}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingQ1<dim, spacedim>::clone() const
{
  return std::make_unique<MappingQ1<dim, spacedim>>(*this);
}

//---------------------------------------------------------------------------


template <int dim, int spacedim>
MappingQGeneric<dim, spacedim>
  StaticMappingQ1<dim, spacedim>::mapping = MappingQGeneric<dim, spacedim>(1);



//--------------------------- Explicit instantiations -----------------------
#include "mapping_q1.inst"


DEAL_II_NAMESPACE_CLOSE
