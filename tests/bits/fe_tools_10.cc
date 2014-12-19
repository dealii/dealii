// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


#include "../tests.h"
#include "fe_tools_common.h"
#include <deal.II/lac/sparsity_pattern.h>

// check
//   FETools::get_fe_from_name
// like fe_tools_09, but with the difference that we use the generic
// "<dim>" marker in the finite element name, instead of the one with
// the concrete dimension (see the documentation)


std::string output_file_name = "output";

template <int dim>
std::string modify_name (const std::string &name)
{
  std::string new_name = name;
  std::string dim_name = std::string("<");
  const char dim_char='0'+dim;
  dim_name += dim_char;
  dim_name += '>';

  std::string::size_type pos;
  while ((pos = new_name.find(dim_name)) != std::string::npos)
    new_name.replace (pos, 3, "<dim>");

  return new_name;
}



template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  FiniteElement<dim> *p1, *p2;

  // check that the name of the fe
  // and the name of the fe that we
  // re-create from this name are
  // identitical. this is also a
  // pretty good indication that the
  // two FEs are actually the same
  deallog << modify_name<dim> (fe1.get_name());
  p1 = FETools::get_fe_from_name<dim> (modify_name<dim> (fe1.get_name()));
  Assert (fe1.get_name() == p1->get_name(),
          ExcInternalError());
  deallog << " ok" << std::endl;
  delete p1;

  // same for fe2
  deallog << modify_name<dim> (fe2.get_name());
  p2 = FETools::get_fe_from_name<dim> (modify_name<dim> (fe2.get_name()));
  Assert (fe2.get_name() == p2->get_name(),
          ExcInternalError());
  deallog << " ok" << std::endl;
  delete p2;
}

