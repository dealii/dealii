// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

// check that we can move a FESystem<dim> object in a reasonable way.

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <string>

#include "../tests.h"

template <int dim>
void
check(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
  deallog << "components: " << fe.n_components() << std::endl;
  deallog << "blocks: " << fe.n_blocks() << std::endl;
  deallog << "conforms H1: " << fe.conforms(FiniteElementData<dim>::H1)
          << std::endl;
  deallog << "n_base_elements: " << fe.n_base_elements() << std::endl;
  deallog << std::endl;
}

template <int dim>
void
move()
{
  FESystem<dim> fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
  check(fe);
  FESystem<dim> fe2(std::move(fe));
  check(fe);
  check(fe2);
}

int
main()
{
  initlog();

  move<1>();
  move<2>();
  move<3>();

  return 0;
}
