// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
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



char logname[] = "output";

#include <deal.II/base/quadrature_lib.h>

#include "injection_common.h"


template <int dim>
void test ()
{
  deallog << std::setprecision (8);

  for (unsigned int i=0; i<4; ++i)
    for (unsigned int j=i; j<4; ++j)
      if (i > 0 && j > 0)
        do_check (FE_DGQArbitraryNodes<dim>(QIterated<1>(QTrapez<1>(),i)),
                  FE_DGQArbitraryNodes<dim>(QIterated<1>(QTrapez<1>(),j)));
      else if (j > 0)
        do_check (FE_DGQ<dim>(0),
                  FE_DGQArbitraryNodes<dim>(QIterated<1>(QTrapez<1>(),j)));
      else
        do_check (FE_DGQ<dim>(0),
                  FE_DGQ<dim>(0));
}
