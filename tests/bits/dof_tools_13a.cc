//----------------------------  dof_tools_13a.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_13a.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.h"
#include <deal.II/lac/vector.h>

// check
//   DoFTools::distribute_cell_to_dof_vector

// this is a variant of dof_tools_13. it turned out that the function
// tested here didn't zero out its argument up front, instead adding
// to the previous content. this is not what we intended, so test with
// a preset vector and make sure that we get the same result as in
// dof_tools_13.



std::string output_file_name = "dof_tools_13a/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
                                   // this doesn't make much sense if
                                   // the element is not primitive
  if (dof_handler.get_fe().is_primitive() == false)
    return;
  
  Vector<double> cell_data (dof_handler.get_tria().n_active_cells());
  for (unsigned int i=0; i<cell_data.size(); ++i)
    cell_data(i) = i;

                                   // distribute to first component
  Vector<double> dof_data (dof_handler.n_dofs());

				   // preset the vector to make sure
				   // that the function zeroes out
				   // previous content. however, only
				   // touch those elements that we
				   // will actually use
  std::vector<bool> component_dofs (dof_handler.n_dofs());
  {
    std::vector<bool> component_mask (dof_handler.get_fe().n_components(),
				      false);
    component_mask[0] = true;
    DoFTools::extract_dofs (dof_handler, component_mask, component_dofs);

    for (unsigned int i=0; i<dof_data.size(); ++i)
      if (component_dofs[i] == true)
	dof_data(i) = i+1;
      else
	dof_data(i) = 0;
  }
  
  DoFTools::distribute_cell_to_dof_vector (dof_handler,
                                           cell_data,
                                           dof_data);
                                   // output every third element
  for (unsigned int i=0; i<dof_data.size(); i+=3)
    deallog << dof_data(i) << " ";
  deallog << std::endl;  

				   // check that no other values were
				   // set
  for (unsigned int i=0; i<dof_data.size(); ++i)
    if (component_dofs[i] == false)
      Assert (dof_data(i) == 0,
	      ExcInternalError());
  

                                   // distribute to last component. by
                                   // default we distribute to
                                   // component zero
  
				   // preset the vector again to make
				   // sure that the function zeroes out
				   // previous content. 
  {
    std::vector<bool> component_mask (dof_handler.get_fe().n_components(),
				      false);
    component_mask.back() = true;
    DoFTools::extract_dofs (dof_handler, component_mask, component_dofs);
    for (unsigned int i=0; i<dof_data.size(); ++i)
      if (component_dofs[i] == true)
	dof_data(i) = i+1;
      else
	dof_data(i) = 0;
  }
  DoFTools::distribute_cell_to_dof_vector (dof_handler,
                                           cell_data,
                                           dof_data,
                                           dof_handler.get_fe().n_components()-1);
  for (unsigned int i=0; i<dof_data.size(); i+=3)
    deallog << dof_data(i) << " ";
  deallog << std::endl;  

				   // check that no other values were
				   // set
  for (unsigned int i=0; i<dof_data.size(); ++i)
    if (component_dofs[i] == false)
      Assert (dof_data(i) == 0,
	      ExcInternalError());
}
