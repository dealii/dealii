// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"

#include "fe_tools_common.h"

// check
//   FETools::get_fe_by_name



template <int dim>
void
check_this(const FiniteElement<dim> &fe1, const FiniteElement<dim> &fe2)
{
  // check that the name of the fe
  // and the name of the fe that we
  // re-create from this name are
  // identitical. this is also a
  // pretty good indication that the
  // two FEs are actually the same
  deallog << fe1.get_name();
  std::unique_ptr<FiniteElement<dim>> p1 =
    FETools::get_fe_by_name<dim, dim>(fe1.get_name());
  AssertThrow(fe1.get_name() == p1->get_name(), ExcInternalError());
  deallog << " ok" << std::endl;

  // same for fe2
  deallog << fe2.get_name();
  std::unique_ptr<FiniteElement<dim>> p2 =
    FETools::get_fe_by_name<dim, dim>(fe2.get_name());
  AssertThrow(fe2.get_name() == p2->get_name(), ExcInternalError());
  deallog << " ok" << std::endl;
}
