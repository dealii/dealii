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
#include "dof_tools_common.h"

// comparison of the results of
//   FE::face_to_equivalent_cell_index (i)
// and
// FE::face_to_cell_index (i, 0)
// According to the documentation they should be the same!


std::string output_file_name = "fe_face_to_cell_index/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  
  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
    {
      deallog << fe.face_to_equivalent_cell_index (i)
	      << ' '
	      << fe.face_to_cell_index (i, 0)
	      << std::endl;
      Assert (fe.face_to_equivalent_cell_index (i) ==
	      fe.face_to_cell_index (i, 0),
	      ExcInternalError());
    }
  
  deallog << "OK" << std::endl;
}

