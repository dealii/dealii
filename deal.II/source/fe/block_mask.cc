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

#include <deal.II/fe/block_mask.h>

#include <iostream>


DEAL_II_NAMESPACE_OPEN

std::ostream & operator << (std::ostream &out,
                            const BlockMask &mask)
{
  if (mask.block_mask.size() == 0)
    out << "[all blocks selected]";
  else
    {
      out << '[';
      for (unsigned int i=0; i<mask.block_mask.size(); ++i)
        {
          out << (mask.block_mask[i] ? "true" : "false");
          if (i != mask.block_mask.size()-1)
            out << ',';
        }
      out << ']';
    }

  return out;
}



std::size_t
BlockMask::memory_consumption () const
{
  return sizeof(*this) + MemoryConsumption::memory_consumption (block_mask);
}


DEAL_II_NAMESPACE_CLOSE
