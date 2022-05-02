// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Check that VectorTools::project_boundary_values_curl_conforming_l2
// also works on systems of Nedelec elements. This is a variation of a
// testcase by Abbas Ballout posted to the mailing list.

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues()
    : Function<dim>(2 * dim)
  {}

  void
  vector_value_list(const std::vector<Point<dim>> &points,
                    std::vector<Vector<double>> &  values) const override
  {
    for (unsigned int i = 0; i < points.size(); ++i)
      {
        values[i][0] = 1.0;
        values[i][1] = 0.0;
        if (dim == 3)
          values[i][2] = 0.0;
        values[i][dim + 0] = 1.0;
        values[i][dim + 1] = 0.0;
        if (dim == 3)
          values[i][dim + 2] = 0.0;
      }
  }
};



template <int dim>
void
test()
{
  deallog << "*** dim=" << dim << std::endl;

  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube(triangulation);
  triangulation.begin_active()
    ->face(dim == 2 ? 2 : 4)
    ->set_boundary_id(1); // set the 4th boundary id to 1 to apply dirichlet BC

  DoFHandler<dim> dof_handler(triangulation);
  FESystem<dim>   fe(FE_Nedelec<dim>(0), 2);
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::component_wise(dof_handler);

  // First try to apply constraints only to the first 3 components
  {
    deallog << "Check 1:" << std::endl;

    AffineConstraints<double> constraints;

    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      0, // starting compenent
      BoundaryValues<dim>(),
      1, // face ID1
      constraints,
      StaticMappingQ1<dim>::mapping);
    constraints.close();
    constraints.print(deallog.get_file_stream());
  }

  // Try again with just the second 3 components
  {
    deallog << "Check 2:" << std::endl;

    AffineConstraints<double> constraints;

    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      dim, // starting component
      BoundaryValues<dim>(),
      1, // face ID1
      constraints,
      StaticMappingQ1<dim>::mapping);
    constraints.close();
    constraints.print(deallog.get_file_stream());
  }

  // And now with all
  {
    deallog << "Check 3:" << std::endl;

    AffineConstraints<double> constraints;

    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      0, // starting compenent
      BoundaryValues<dim>(),
      1, // face ID1
      constraints,
      StaticMappingQ1<dim>::mapping);
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      dim, // starting compenent
      BoundaryValues<dim>(),
      1, // face ID1
      constraints,
      StaticMappingQ1<dim>::mapping);
    constraints.close();
    constraints.print(deallog.get_file_stream());
  }
}



int
main()
{
  initlog();

  test<2>();
  test<3>();
}
