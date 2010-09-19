//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// like the _01 test, except that we declare output component fields
// as vectors if possible. since VTK is the only format that currently
// supports this, use it as the only output format. (deal.II
// intermediate also supports vector-valued components, but we already
// test for this in the _04 tests)
//
// attach the same data twice to see if we can produce output with
// more than one vector field

#include "../tests.h"
#include "data_out_common.h"
#include <lac/sparsity_pattern.h>
#include <numerics/data_out.h>


std::string output_file_name = "data_out_05/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler,
            const Vector<double>  &v_node,
            const Vector<double>  &v_cell)
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation (dof_handler.get_fe().n_components(),
				   DataComponentInterpretation::component_is_scalar);
				   // if possible, declare the last
				   // dim components as vectors
  if (dof_handler.get_fe().n_components() >= dim)
    for (unsigned int i=dof_handler.get_fe().n_components()-dim;
	 i<dof_handler.get_fe().n_components(); ++i)
      data_component_interpretation[i]
	= DataComponentInterpretation::component_is_part_of_vector;

  deallog << "Checking " << dof_handler.get_fe().get_name()
	  << std::endl
	  << "Component mask: ";
  for (unsigned int i=0; i<dof_handler.get_fe().n_components(); ++i)
    deallog << (data_component_interpretation[i] ==
		DataComponentInterpretation::component_is_scalar
		?
		"scalar "
		:
		"vector ");
  deallog << std::endl;
  
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (v_node,
			    std::vector<std::string>(dof_handler.get_fe().n_components(),
						     "node_data_1"),
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  data_out.add_data_vector (v_node,
			    std::vector<std::string>(dof_handler.get_fe().n_components(),
						     "node_data_2"),
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  data_out.add_data_vector (v_cell, "cell_data");
  data_out.build_patches ();
  
  data_out.write_vtk (deallog.get_file_stream());
}


