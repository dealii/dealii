// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
#include "hp/mapping_collection.inst"


DEAL_II_NAMESPACE_CLOSE
