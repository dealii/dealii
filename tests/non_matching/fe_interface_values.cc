// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/non_matching/fe_values.h>

#include <deal.II/numerics/vector_tools.h>

#include <optional>

#include "../tests.h"


// Set up a triangulation with 2 elements in a row: |-0-|-1-|, over
// [0,2]x[0,1]^{dim-1}, and a level set  function with a cut in the plane
// $x_[dim-1] = 0.5$,  such that it cuts through the middle of both cells in 2D
// and 3D and thought the middle of cell 0 in 1D. Test that we can create an
// object of type NonMatching::FEInterfaceValues. Call reinit on all the faces
// of cell 0 and check that the dealii::FEInterfaceValues objects we get from
// get_inside/outside_fe_values are set up correctly. That is, the
// std::optionals that should contain a value do and the ones that should not
// contain a value do not.
template <int dim>
class Test
{
public:
  Test();

  void
  run();

private:
  void
  setup_mesh();

  void
  setup_discrete_level_set();

  template <typename IteratorType>
  void
  print_which_optionals_have_values_on_cell_0(
    NonMatching::FEInterfaceValues<dim> &fe_values,
    IteratorType                         cell);

  Triangulation<dim>    triangulation;
  hp::FECollection<dim> fe_collection;
  DoFHandler<dim>       dof_handler;

  hp::MappingCollection<dim>           mapping_collection;
  hp::QCollection<dim - 1>             q_collection;
  hp::QCollection<1>                   q_collection1D;
  const NonMatching::RegionUpdateFlags region_update_flags;

  Vector<double>                   level_set;
  NonMatching::MeshClassifier<dim> mesh_classifier;
};



template <int dim>
Test<dim>::Test()
  : dof_handler(triangulation)
  , mesh_classifier(dof_handler, level_set)
{
  fe_collection.push_back(FE_Q<dim>(1));
  mapping_collection.push_back(MappingCartesian<dim>());

  const unsigned int n_quadrature_points = 1;
  q_collection.push_back(QGauss<dim - 1>(n_quadrature_points));
  q_collection1D.push_back(QGauss<1>(n_quadrature_points));
}



template <int dim>
void
Test<dim>::run()
{
  setup_mesh();
  setup_discrete_level_set();
  mesh_classifier.reclassify();

  // Test to create an object with each constructor.
  {
    deallog << "Simple constructor:" << std::endl;
    NonMatching::FEInterfaceValues<dim> fe_values(fe_collection,
                                                  q_collection1D[0],
                                                  region_update_flags,
                                                  mesh_classifier,
                                                  dof_handler,
                                                  level_set);
    print_which_optionals_have_values_on_cell_0(fe_values,
                                                triangulation.begin_active());
    print_which_optionals_have_values_on_cell_0(fe_values,
                                                dof_handler.begin_active());
  }
  {
    deallog << "Advanced constructor:" << std::endl;
    NonMatching::FEInterfaceValues<dim> fe_values(mapping_collection,
                                                  fe_collection,
                                                  q_collection,
                                                  q_collection1D,
                                                  region_update_flags,
                                                  mesh_classifier,
                                                  dof_handler,
                                                  level_set);
    print_which_optionals_have_values_on_cell_0(fe_values,
                                                triangulation.begin_active());
    print_which_optionals_have_values_on_cell_0(fe_values,
                                                dof_handler.begin_active());
  }
}



template <int dim>
void
Test<dim>::setup_mesh()
{
  const Point<dim> lower_left;
  Point<dim>       upper_right;

  std::vector<unsigned int> repetitions;
  upper_right[0] = 2;
  repetitions.push_back(2);
  for (unsigned int d = 1; d < dim; ++d)
    {
      upper_right[d] = 1;
      repetitions.push_back(1);
    }

  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            repetitions,
                                            lower_left,
                                            upper_right);
}



template <int dim>
void
Test<dim>::setup_discrete_level_set()
{
  dof_handler.distribute_dofs(fe_collection);
  level_set.reinit(dof_handler.n_dofs());

  Point<dim> point_on_zero_contour;
  point_on_zero_contour[dim - 1] = 0.5;
  const Functions::SignedDistance::Plane<dim> analytical_levelset(
    point_on_zero_contour, Point<dim>::unit_vector(dim - 1));

  VectorTools::interpolate(dof_handler, analytical_levelset, level_set);
}



template <int dim>
template <typename IteratorType>
void
Test<dim>::print_which_optionals_have_values_on_cell_0(
  NonMatching::FEInterfaceValues<dim> &fe_values,
  IteratorType                         cell)
{
  for (const unsigned int face_index : cell->face_indices())
    {
      deallog << "face " << face_index << std::endl;

      if (cell->at_boundary(face_index))
        {
          fe_values.reinit(cell, face_index);
        }
      else
        {
          const unsigned int invalid_subface =
            dealii::numbers::invalid_unsigned_int;

          fe_values.reinit(cell,
                           face_index,
                           invalid_subface,
                           cell->neighbor(face_index),
                           cell->neighbor_of_neighbor(face_index),
                           invalid_subface);
        }

      if (fe_values.get_inside_fe_values())
        deallog << "inside has value" << std::endl;

      if (fe_values.get_outside_fe_values())
        deallog << "outside has value" << std::endl;
    }
}



template <int dim>
void
run_test()
{
  deallog << "dim = " << dim << std::endl;
  Test<dim> test;
  test.run();
  deallog << std::endl;
}



int
main()
{
  initlog();

  run_test<1>();
  run_test<2>();
  run_test<3>();
}
