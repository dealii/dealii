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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// check that we can move a FE_Q<dim> object in a reasonable way.

#include <deal.II/fe/fe_q.h>

#include <string>

#include "../tests.h"

template <int dim>
void
check(const FiniteElement<dim> &fe)
{
  deallog << "dim: " << dim << std::endl;
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
  FE_Q<dim> fe(1);
  check(fe);
  FE_Q<dim> fe2(std::move(fe));
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
