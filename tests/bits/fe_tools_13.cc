// ---------------------------------------------------------------------
// $Id$
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
//   FE::face_to_equivalent_cell_index


std::string output_file_name = "fe_tools_13/output";



template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  for (unsigned int i=0; i<fe1.dofs_per_face; ++i)
    deallog << fe1.face_to_equivalent_cell_index(i) << std::endl;

  for (unsigned int i=0; i<fe2.dofs_per_face; ++i)
    deallog << fe2.face_to_equivalent_cell_index(i) << std::endl;
}

