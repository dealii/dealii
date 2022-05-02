// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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



// Loop over entries of FECollection.

#include <deal.II/fe/fe_q.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int dim = 2;

  hp::FECollection<dim> collection(FE_Q<dim>(1), FE_Q<dim>(2), FE_Q<dim>(3));

  for (const auto &fe : collection)
    deallog << fe.n_dofs_per_cell() << std::endl;

  return 0;
}
