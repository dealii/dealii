// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
