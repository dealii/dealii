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



// this is reduced from hp_constraints_q_system_x_01, where we test that we
// can deal with FESystem(FE_Q(p),FE_DGQ(q)) for different p,q. note that for
// fixed p but varying q, neither of the two elements will dominate the other

char logname[] = "output";


#include "hp_constraints_common.h"


template <int dim>
void test ()
{
  hp::FECollection<dim> fe;
  fe.push_back (FESystem<dim>(FE_Q<dim>(1), 1,
                              FE_DGQ<dim>(0), 1));
  fe.push_back (FESystem<dim>(FE_Q<dim>(1), 1,
                              FE_DGQ<dim>(1), 1));

  test_no_hanging_nodes (fe);
}
