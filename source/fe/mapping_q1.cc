// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/mapping_q1.h>


DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
MappingQ1<dim, spacedim>::MappingQ1()
  : MappingQ<dim, spacedim>(1)
{}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingQ1<dim, spacedim>::clone() const
{
  return std::make_unique<MappingQ1<dim, spacedim>>(*this);
}

//--------------------------- Explicit instantiations -----------------------
#include "fe/mapping_q1.inst"


DEAL_II_NAMESPACE_CLOSE
