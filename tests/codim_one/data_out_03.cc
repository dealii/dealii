// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// DataOut::build_patches did not compute mapped coordinates for all cells
// if dim<spacedim even if a mapping was explicitly requested.

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
class Identity : public Function<dim>
{
public:
  Identity()
    : Function<dim>(dim)
  {}


  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return p[component];
  }

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const
  {
    for (unsigned int i = 0; i < dim; ++i)
      values(i) = p[i];
  }
};


int
main()
{
  initlog();

  int fe_degree      = 2;
  int mapping_degree = 2;


  Triangulation<2, 3> tria;

  std::map<Triangulation<2, 3>::cell_iterator,
           Triangulation<3, 3>::face_iterator>
    surface_to_volume_mapping;

  SphericalManifold<3> boundary_description;
  Triangulation<3>     volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);


  volume_mesh.set_manifold(1, boundary_description);
  volume_mesh.set_manifold(0, boundary_description);

  static SphericalManifold<2, 3> surface_description;
  tria.set_manifold(1, surface_description);
  tria.set_manifold(0, surface_description);

  const std::set<types::boundary_id> boundary_ids = {0};
  GridGenerator::extract_boundary_mesh(volume_mesh, tria, boundary_ids);

  // test for the position
  MappingQ<2, 3> mapping(mapping_degree);

  FESystem<2, 3>   fe_test(FE_Q<2, 3>(fe_degree), 3);
  DoFHandler<2, 3> dh_test(tria);
  dh_test.distribute_dofs(fe_test);

  Vector<double> position(dh_test.n_dofs());
  VectorTools::interpolate(mapping, dh_test, Identity<3>(), position);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      3, DataComponentInterpretation::component_is_part_of_vector);

  std::vector<std::string> solution_names(3, "position");

  DataOut<2, 3> data_out;
  data_out.attach_dof_handler(dh_test);
  data_out.add_data_vector(position,
                           solution_names,
                           DataOut<2, 3>::type_dof_data,
                           data_component_interpretation);
  data_out.build_patches(mapping, 2);

  data_out.write_gnuplot(deallog.get_file_stream());


  dh_test.clear();
  tria.reset_manifold(0);

  return 0;
}
