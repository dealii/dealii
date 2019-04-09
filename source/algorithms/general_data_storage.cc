// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


#include <deal.II/algorithms/general_data_storage.h>


DEAL_II_NAMESPACE_OPEN


std::size_t
GeneralDataStorage::size() const
{
  return any_data.size();
}


void
GeneralDataStorage::merge(const GeneralDataStorage &other)
{
  any_data.insert(other.any_data.begin(), other.any_data.end());
}


void
GeneralDataStorage::reset()
{
  any_data.clear();
}


bool
GeneralDataStorage::stores_object_with_name(const std::string &name) const
{
  return any_data.find(name) != any_data.end();
}


void
GeneralDataStorage::remove_object_with_name(const std::string &name)
{
  const auto it = any_data.find(name);
  if (it != any_data.end())
    any_data.erase(it);
}



DEAL_II_NAMESPACE_CLOSE
