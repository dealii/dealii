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

// Check that FiniteElement::get_generalized_support_points() returns a
// vector of unique generalized support points.

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <vector>

#include "../tests.h"


template <int dim>
void
test()
{
  const auto print = [](const FiniteElement<dim> &fe) {
    deallog << fe.get_name() << std::endl;
    for (auto it : fe.get_generalized_support_points())
      {
        deallog << "(" << it[0] << "," << it[1] << ") ";
      }
    deallog << std::endl;
  };

  deallog << "dim " << dim << " lowest order" << std::endl;
  print(FE_Q<dim>(1));
  print(FE_RaviartThomas<dim>(0));
  print(FE_Nedelec<dim>(0));
  print(FESystem<dim>(FE_RaviartThomas<dim>(0), 3, FE_Nedelec<dim>(0), 3));
  print(FESystem<dim>(
    FESystem<dim>(FE_RaviartThomas<dim>(0), 3, FE_Nedelec<dim>(0), 3),
    2,
    FE_Q<dim>(1),
    3,
    FE_Q<dim>(1),
    3));

  deallog << "dim " << dim << " higher order" << std::endl;
  print(FE_Q<dim>(2));
  print(FE_RaviartThomas<dim>(1));
  print(FE_Nedelec<dim>(1));
  print(FESystem<dim>(FE_RaviartThomas<dim>(1), 3, FE_Nedelec<dim>(1), 3));
  print(FESystem<dim>(
    FESystem<dim>(FE_RaviartThomas<dim>(1), 3, FE_Nedelec<dim>(1), 3),
    2,
    FE_Q<dim>(2),
    3,
    FE_Q<dim>(1),
    3));
}

int
main()
{
  initlog();

  test<2>();
  test<3>();
}
