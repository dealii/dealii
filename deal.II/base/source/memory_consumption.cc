//----------------------------  function.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function.cc  ---------------------------


#include <base/memory_consumption.h>


namespace MemoryConsumption 
{
  unsigned int memory_consumption (const std::vector<std::string> &v)
  {
    unsigned int mem = sizeof(v);
    for (unsigned int i=0; i<v.size(); ++i)
      mem += memory_consumption(v[i]);
    return mem;
  };
  

};
