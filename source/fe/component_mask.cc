// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/fe/component_mask.h>

#include <iostream>


DEAL_II_NAMESPACE_OPEN

std::ostream &operator << (std::ostream &out,
                           const ComponentMask &mask)
{
  if (mask.component_mask.size() == 0)
    out << "[all components selected]";
  else
    {
      out << '[';
      for (unsigned int i=0; i<mask.component_mask.size(); ++i)
        {
          out << (mask.component_mask[i] ? "true" : "false");
          if (i != mask.component_mask.size()-1)
            out << ',';
        }
      out << ']';
    }

  return out;
}



std::size_t
ComponentMask::memory_consumption () const
{
  return sizeof(*this) + MemoryConsumption::memory_consumption (component_mask);
}


DEAL_II_NAMESPACE_CLOSE
