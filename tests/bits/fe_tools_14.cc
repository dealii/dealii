// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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

// check
//   FE::hp_constraints_are_implemented ()
// a bit like hp_constraints_are_implemented, but with a different
// set of elements


std::string output_file_name = "output";



template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  deallog << (fe1.hp_constraints_are_implemented () ? "true" : "false")
          << std::endl;
  deallog << (fe2.hp_constraints_are_implemented () ? "true" : "false")
          << std::endl;
}

