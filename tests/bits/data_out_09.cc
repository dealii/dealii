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


// check that add_data_vector and type_automatic work correctly and confirm
// that the Assert is triggered that checks, when we have dof_data and
// cell_data with the same length.


#include "../tests.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>




void test()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global(1);
  
  FE_DGQ<2> fe(0);
  DoFHandler<2> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  std::vector<types::global_dof_index> renumbering(dof_handler.n_dofs());
  renumbering[0]=3;
  renumbering[1]=2;
  renumbering[2]=1;
  renumbering[3]=0;
  
  dof_handler.renumber_dofs(renumbering);
  

  Vector<double> v_node(dof_handler.n_dofs());
  Vector<double> v_cell(dof_handler.n_dofs());
  for (unsigned int i=0;i<dof_handler.n_dofs();++i)
    {
      v_node[i]=i;
      v_cell[i]=10+i;
    }

  //correct one is done if both are possible but specified
  deallog << "*** Check cell and node data:" << std::endl;
  {    
    DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (v_node, "node_data", DataOut<2>::type_dof_data);
    data_out.add_data_vector (v_cell, "cell_data", DataOut<2>::type_cell_data);
    data_out.build_patches ();

    data_out.write_vtk (deallog.get_file_stream());
  }

  //only tria, output correctly
  deallog << "*** Check cell data only:" << std::endl;
  {  
    DataOut<2> data_out;
    data_out.attach_triangulation (tria);
    data_out.add_data_vector (v_cell, "cell_data");
    data_out.build_patches ();

    data_out.write_vtk (deallog.get_file_stream());
  }
  

  //error if both
  deallog << "*** should fail:" << std::endl;
  try
  {  
    DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (v_cell, "cell_data");
    data_out.build_patches ();

    data_out.write_vtk (deallog.get_file_stream());
  }
  catch (...)
    {
      deallog << "exception" << std::endl;
    }

  // works if we use the other add_data_vector:
  deallog << "*** Check node data only:" << std::endl;
  {  
    DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (dof_handler, v_node, "node_data");
    data_out.build_patches ();

    data_out.write_vtk (deallog.get_file_stream());
  }
  
  
  
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();
  initlog();
  test();
  
}

