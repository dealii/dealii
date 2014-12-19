// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


/*
 * Project the function [1,1] onto a deformed grid and see whether the ABF
 * elements can represent it exactly.
 */



char logname[] = "output";
#include "deformed_projection.h"


void test ()
{
  FE_ABF<2> fe (0);
  check (fe);
}
