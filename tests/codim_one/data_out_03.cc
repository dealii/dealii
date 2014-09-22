// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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



// DataOut::build_patches did not compute mapped coordinates for all cells
// if dim<spacedim even if a mapping was explicitly requested.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/mapping_q.h>

std::ofstream logfile("output");

template <int dim>
class Identity : public Function<dim>
{
public:
  Identity ()
    :
    Function<dim>(dim)
  {
  }


  virtual double value (const Point<dim> &p,
                        const unsigned int component) const
  {
    return p(component);
  }

  virtual void vector_value(const Point<dim> &p,
                            Vector< double > &values) const
  {
    for (unsigned int i=0; i<dim; i++)
      values(i) = p(i);
  }
};


int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  int fe_degree =2;
  int mapping_degree = 2;


  Triangulation<2,3> tria;

  std::map< Triangulation<2,3>::cell_iterator,
      Triangulation<3,3>::face_iterator> surface_to_volume_mapping;

  HyperBallBoundary<3> boundary_description;
  Triangulation<3> volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);


  volume_mesh.set_boundary (1, boundary_description);
  volume_mesh.set_boundary (0, boundary_description);

  static HyperBallBoundary<2,3> surface_description;
  tria.set_boundary (1, surface_description);
  tria.set_boundary (0, surface_description);

  std::set<types::boundary_id> boundary_ids;
  boundary_ids.insert(0);

  GridGenerator::extract_boundary_mesh (volume_mesh, tria,
                                    boundary_ids);

  // test for the position
  MappingQ<2,3>   mapping(mapping_degree,true);

  FESystem<2,3> fe_test (FE_Q<2,3>
                         (fe_degree),3);
  DoFHandler<2,3> dh_test(tria);
  dh_test.distribute_dofs(fe_test);

  Vector<double> position(dh_test.n_dofs());
  VectorTools::interpolate(mapping,dh_test,Identity<3>(),position);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  data_component_interpretation
  (3, DataComponentInterpretation::component_is_part_of_vector);

  std::vector<std::string> solution_names (3, "position");

  DataOut<2, DoFHandler<2,3> > data_out;
  data_out.attach_dof_handler (dh_test);
  data_out.add_data_vector (position, solution_names,
                            DataOut<2,DoFHandler<2,
                            3> >::type_dof_data,
                            data_component_interpretation);
  data_out.build_patches (mapping, 2);

  data_out.write_gnuplot (deallog.get_file_stream());


  dh_test.clear();
  tria.set_boundary(0);

  return 0;
}
