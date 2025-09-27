// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"

#include "fe_tools_common.h"

// check
//   FETools::get_fe_by_name
// like fe_tools_09 and fe_tools_10, but this time with no dimension
// marker at all (see the documentation)



template <int dim>
std::string
modify_name(const std::string &name)
{
  std::string new_name = name;
  std::string dim_name = std::string("<");
  const char  dim_char = '0' + dim;
  dim_name += dim_char;
  dim_name += '>';

  std::string::size_type pos;
  while ((pos = new_name.find(dim_name)) != std::string::npos)
    new_name.replace(pos, 3, "");

  return new_name;
}



template <int dim>
void
check_this(const FiniteElement<dim> &fe1, const FiniteElement<dim> &fe2)
{
  // check that the name of the FE
  // and the name of the FE that we
  // re-create from this name are
  // identical. this is also a
  // pretty good indication that the
  // two FEs are actually the same
  deallog << modify_name<dim>(fe1.get_name());
  std::unique_ptr<FiniteElement<dim>> p1 =
    FETools::get_fe_by_name<dim, dim>(modify_name<dim>(fe1.get_name()));
  AssertThrow(fe1.get_name() == p1->get_name(), ExcInternalError());
  deallog << " ok" << std::endl;

  // same for fe2
  deallog << modify_name<dim>(fe2.get_name());
  std::unique_ptr<FiniteElement<dim>> p2 =
    FETools::get_fe_by_name<dim, dim>(modify_name<dim>(fe2.get_name()));
  AssertThrow(fe2.get_name() == p2->get_name(), ExcInternalError());
  deallog << " ok" << std::endl;
}
