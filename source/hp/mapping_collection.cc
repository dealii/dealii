// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2021 by the deal.II authors
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


#include <deal.II/base/memory_consumption.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/hp/mapping_collection.h>

DEAL_II_NAMESPACE_OPEN


namespace hp
{
  template <int dim, int spacedim>
  MappingCollection<dim, spacedim>::MappingCollection(
    const Mapping<dim, spacedim> &mapping)
  {
    this->push_back(mapping);
  }



  template <int dim, int spacedim>
  MappingCollection<dim, spacedim>::MappingCollection(
    const MappingCollection<dim, spacedim> &other)
  {
    for (unsigned int i = 0; i < other.size(); ++i)
      push_back(other[i]);
  }



  template <int dim, int spacedim>
  void
  MappingCollection<dim, spacedim>::push_back(
    const Mapping<dim, spacedim> &new_mapping)
  {
    Collection<Mapping<dim, spacedim>>::push_back(
      std::shared_ptr<const Mapping<dim, spacedim>>(new_mapping.clone()));
  }

} // namespace hp



// explicit instantiations
#include "mapping_collection.inst"


DEAL_II_NAMESPACE_CLOSE
