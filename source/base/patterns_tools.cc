// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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


#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns_tools.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/std_cxx14/memory.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/core/demangle.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

//#include <fstream>
//#include <iostream>
//#include <iomanip>
//#include <cstdlib>
//#include <algorithm>
//#include <sstream>
//#include <cctype>
//#include <limits>
//#include <cstring>


DEAL_II_NAMESPACE_OPEN



//TODO[WB]: various functions here could be simplified by using namespace Utilities

namespace Patterns
{


  std::string default_list_separator(unsigned int rank)
  {
    static std::array<std::string, 5> seps = {{" ",  ","  ,  ";"  ,  "|"  ,   "%"}};
    AssertThrow(rank < seps.size(), ExcMessage("I don't know what to use for such "
                                               "high rank. Bailing out."));
    return seps[rank];
  }
}

DEAL_II_NAMESPACE_CLOSE
