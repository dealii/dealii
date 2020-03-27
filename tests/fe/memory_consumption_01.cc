// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Test memory_consumption of FE_Q and FE_DGQ.


#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include "../tests.h"


template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << " " << fe.memory_consumption() << std::endl;
}

template <int dim>
void
tests()
{
  // FE_Q
  for (unsigned int i = 1; i < 4; i++)
    test(FE_Q<dim>(i));

  // FE_DGQ
  for (unsigned int i = 0; i < 4; i++)
    test(FE_DGQ<dim>(i));
}

int
main()
{
  initlog();

  {
    deallog.push("1d");
    tests<1>();
    deallog.pop();
  }
  {
    deallog.push("2d");
    tests<2>();
    deallog.pop();
  }
  {
    deallog.push("3d");
    tests<3>();
    deallog.pop();
  }
}
