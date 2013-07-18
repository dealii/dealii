// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2013 by the deal.II authors
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

// This file compiles the second half of the instantiations from solution_transfer.cc
// to get the memory consumption below 1.5gb with gcc (if compiling with PETSc and Trilinos).

#define SOLUTION_TRANSFER_INSTANTIATE_PART_TWO

#include "solution_transfer.cc"
