// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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


#include <deal.II/base/parallel.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace Vector
  {
    // set minimum grain size. this value is
    // roughly in accordance with the curve
    // in the TBB book (fig 3.2) that shows
    // run time as a function of grain size
    // -- there, values from 200 upward are
    // so that the scheduling overhead
    // amortizes well (for very large values
    // in that example, the grain size is too
    // large to split the work load into
    // enough chunks and the problem becomes
    // badly balanced)
    unsigned int minimum_parallel_grain_size = 1000;
  }


  namespace SparseMatrix
  {
    // set this value to 1/5 of the value of
    // the minimum grain size of
    // vectors. this rests on the fact that
    // we have to do a lot more work per row
    // of a matrix than per element of a
    // vector. it could possibly be reduced
    // even further but that doesn't appear
    // worth it any more for anything but
    // very small matrices that we don't care
    // that much about anyway.
    unsigned int minimum_parallel_grain_size = 200;
  }
}




DEAL_II_NAMESPACE_CLOSE
