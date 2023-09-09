// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
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


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include "../tests.h"


/**
 * Test that we can call get_update_flags() on FEInterfaceValues. Test added
 * because this previously led to an address boundary error when using the hp
 * capabilities.
 */
template <int dim>
void
test()
{
  const MappingQ<dim>              mapping(1);
  const hp::MappingCollection<dim> mapping_collection(mapping);

  const FE_Q<dim>             fe(1);
  const hp::FECollection<dim> fe_collection(fe);

  const QGauss<dim - 1>          quadrature(fe.degree + 1);
  const hp::QCollection<dim - 1> q_collection(quadrature);

  const UpdateFlags update_flags = update_default;

  {
    FEInterfaceValues<dim> fiv(mapping_collection,
                               fe_collection,
                               q_collection,
                               update_flags);
    fiv.get_update_flags();
  }
  {
    FEInterfaceValues<dim> fiv(mapping, fe, quadrature, update_flags);
    fiv.get_update_flags();
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();
  test<2>();
}
