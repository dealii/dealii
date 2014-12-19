// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// check that VectorTools::project works for RT elements correctly

char logname[] = "output";


#include "project_common.h"


template <int dim>
void test ()
{
  if (dim != 1)
    // this is interesting also in 3d, but is
    // exceedingly slow there. limit to the
    // case of RT(0) elements in 3d
    for (unsigned int p=0; p<(dim == 2 ? 3 : 1); ++p)
      test_with_2d_deformed_refined_mesh (FE_RaviartThomas<dim>(p), p+1, 1);
}
