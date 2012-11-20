//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

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
