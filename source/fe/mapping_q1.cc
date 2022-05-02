// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

//---------------------------------------------------------------------------


template <int dim, int spacedim>
MappingQ<dim, spacedim>
  StaticMappingQ1<dim, spacedim>::mapping = MappingQ<dim, spacedim>(1);



//--------------------------- Explicit instantiations -----------------------
#include "mapping_q1.inst"


DEAL_II_NAMESPACE_CLOSE
