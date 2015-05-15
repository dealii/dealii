// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2014 by the deal.II authors
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



// DataOut::build_patches appeared to have a problem with outputting
// lines in 2d where nodes were numbered differently when writing data
// vectors as opposed to writing node locations. in the end this
// turned out to be a feature: the mesh was a circle of lines, so
// there are equally many cells as their were nodes, and consequently
// DataOut assumed that it had cell_data, rather than
// dof_data. passing the correct argument fixed the problem, but it
// won't hurt to have this test anyway.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/fe/mapping_fe_field.h>

std::ofstream logfile("output");


int main ()
{
  initlog();

  const unsigned int dim = 2;
  const unsigned int degree = 3;
  const unsigned int n_cycles = 4;
  Triangulation<dim, dim+1>    triangulation;


  FE_Q<dim, dim+1>             fe(degree);
  FESystem<dim, dim+1>         fe_euler(FE_Q<dim, dim+1> (degree), dim+1);

  DoFHandler<dim, dim+1>       dof_handler(triangulation);
  DoFHandler<dim, dim+1>       map_dh(triangulation);


  Mapping<dim, dim+1>  *mapping;

  Vector<double> euler_vec;
  Vector<double> scal_sol;

  GridIn<dim, dim+1> grid_in;
  grid_in.attach_triangulation (triangulation);
  std::ifstream fname(SOURCE_DIR "/grids/sphere_0.inp");
  grid_in.read_ucd (fname);

  SphericalManifold<dim,dim+1> manifold;
  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, manifold);


  for(unsigned int cycle=0; cycle<n_cycles; ++cycle)
  {
    dof_handler.distribute_dofs (fe);
    map_dh.distribute_dofs (fe_euler);

    euler_vec.reinit (map_dh.n_dofs());
    scal_sol.reinit (dof_handler.n_dofs());
    scal_sol = 1;
    std::cout<<euler_vec.size()<<" "<<map_dh.n_dofs()<<std::endl;
    VectorTools::get_position_vector(map_dh,euler_vec);
    if(cycle == 0)
      mapping = new MappingFEField<dim, dim+1>(map_dh,euler_vec);

    if(cycle != n_cycles-1)
      triangulation.refine_global();

  }



  DataOut<dim, DoFHandler<dim, dim+1> > data_out_scal;
  data_out_scal.attach_dof_handler (dof_handler);

  data_out_scal.add_data_vector (scal_sol, "scalar_data",
                            DataOut<dim,DoFHandler<dim, dim+1> >::type_dof_data);

  data_out_scal.build_patches(*mapping,
                              degree,
                              DataOut<dim, DoFHandler<dim, dim+1> >::curved_inner_cells);

  std::string filename_scal = ( "scal_check_"+ Utilities::int_to_string(degree) +
                     ".vtu" );
  std::ofstream file_scal(filename_scal.c_str());

  data_out_scal.write_vtu(file_scal);

  data_out_scal.write_vtk(deallog.get_file_stream());


  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  data_component_interpretation
  (dim+1, DataComponentInterpretation::component_is_part_of_vector);

  DataOut<dim, DoFHandler<dim, dim+1> > data_out_euler;
  data_out_euler.attach_dof_handler (map_dh);

  data_out_euler.add_data_vector (euler_vec, "euler_vec",
        DataOut<dim, DoFHandler<dim, dim+1> >::type_dof_data, data_component_interpretation);
  data_out_euler.build_patches(*mapping,
                        degree,
                        DataOut<dim, DoFHandler<dim, dim+1> >::curved_inner_cells);

  std::string filename_euler = ( "euler_check_"+ Utilities::int_to_string(degree) +
                     ".vtu" );
  std::ofstream file_euler(filename_euler.c_str());

  data_out_euler.write_vtu(file_euler);

  data_out_euler.write_vtk(deallog.get_file_stream());


  if(mapping)
    delete(mapping);

  triangulation.set_manifold(0);

  return 0;
}
