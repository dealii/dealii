// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// this function tests whether the compression of mapping (Jacobians) works
// properly. There should only be a few different Jacobians also when there
// are many cells as the weights should be identical

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"

#include "create_mesh.h"



template <int dim>
void
test()
{
  deallog << "General mesh" << std::endl;
  Triangulation<dim> tria;
  create_mesh(tria);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  typename Triangulation<dim>::active_cell_iterator cell, endc;
  cell = tria.begin_active();
  endc = tria.end();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 0.5)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  for (unsigned int i = 0; i < 5 - dim; ++i)
    {
      cell                 = tria.begin_active();
      endc                 = tria.end();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (cell->center()[0] < 5.)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  const QGauss<1>                          quad(2);
  MatrixFree<dim>                          mf;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  mf.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);

  const unsigned int        n_cell_batches = mf.n_cell_batches();
  std::vector<unsigned int> n_cell_types(4, 0);
  for (unsigned int i = 0; i < n_cell_batches; ++i)
    n_cell_types[mf.get_mapping_info().get_cell_type(i)]++;

  // should do at least some compression
  Assert(n_cell_types[0] + n_cell_types[1] > 0, ExcInternalError());
  Assert(mf.get_mapping_info().cell_data[0].jacobians[0].size() <
           (n_cell_types[3] * quad.size() + n_cell_batches - n_cell_types[3]),
         ExcInternalError());
  deallog << "OK" << std::endl;
}



template <int dim>
void
test_cube()
{
  deallog << "Hyper cube" << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5 - dim);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  const QGauss<1>                          quad(2);
  MatrixFree<dim>                          mf;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  mf.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);

  const unsigned int        n_cell_batches = mf.n_cell_batches();
  std::vector<unsigned int> n_cell_types(4, 0);
  for (unsigned int i = 0; i < n_cell_batches; ++i)
    n_cell_types[mf.get_mapping_info().get_cell_type(i)]++;

  // should have one Cartesian cell and no other cell type
  AssertDimension(n_cell_types[0], n_cell_batches);
  AssertDimension(mf.get_mapping_info().cell_data[0].jacobians[0].size(), 2);
  Assert(n_cell_batches > 1, ExcInternalError());
  deallog << "OK" << std::endl;
}



template <int dim>
void
test_parallelogram()
{
  deallog << "Parallelogram" << std::endl;
  Triangulation<dim> tria;
  Point<dim>         corners[dim];
  for (unsigned int d = 0; d < dim; ++d)
    {
      corners[d][d] = 1.;
      if (d > 0)
        corners[d][0] = 0.5 + 0.5 * d;
    }
  GridGenerator::parallelepiped(tria, corners);
  tria.refine_global(5 - dim);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  const QGauss<1>                          quad(2);
  MatrixFree<dim>                          mf;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  mf.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);

  const unsigned int        n_cell_batches = mf.n_cell_batches();
  std::vector<unsigned int> n_cell_types(4, 0);
  for (unsigned int i = 0; i < n_cell_batches; ++i)
    n_cell_types[mf.get_mapping_info().get_cell_type(i)]++;

  // should have one affine cell and no other
  // cell type
  AssertDimension(n_cell_types[1], n_cell_batches);
  AssertDimension(mf.get_mapping_info().cell_data[0].jacobians[0].size(), 2);
  Assert(n_cell_batches > 1, ExcInternalError());
  deallog << "OK" << std::endl;
}



template <int dim>
class DeformedManifold : public ChartManifold<dim, dim, dim>
{
public:
  DeformedManifold() = default;

  virtual std::unique_ptr<Manifold<dim, dim>>
  clone() const
  {
    return std::make_unique<DeformedManifold<dim>>();
  }

  virtual Point<dim>
  push_forward(const Point<dim> &chart_point) const
  {
    Point<dim> result = chart_point;
    result[0]         = std::tanh(2.0 * result[0]) / std::tanh(2.0);
    return result;
  }

  virtual Point<dim>
  pull_back(const Point<dim> &space_point) const
  {
    Point<dim> result = space_point;
    result[0]         = 0.5 * std::atanh(result[0] * std::tanh(2.0));
    return result;
  }
};



template <int dim>
void
test_deformed_cube()
{
  deallog << "Deformed hyper cube" << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1., 1.);
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, DeformedManifold<dim>());
  tria.refine_global(6 - dim);

  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  constraints.close();

  const QGauss<1>                          quad(3);
  MatrixFree<dim>                          mf;
  typename MatrixFree<dim>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim>::AdditionalData::none;
  data.mapping_update_flags_inner_faces =
    update_gradients | update_normal_vectors;
  data.mapping_update_flags_boundary_faces =
    update_gradients | update_normal_vectors;

  mf.reinit(MappingQ1<dim>{}, dof, constraints, quad, data);
  const unsigned int n_cell_batches = mf.n_cell_batches();
  Assert(n_cell_batches > 1, ExcInternalError());

  {
    std::vector<unsigned int> n_cell_types(4, 0);
    for (unsigned int i = 0; i < n_cell_batches; ++i)
      n_cell_types[mf.get_mapping_info().get_cell_type(i)]++;

    // should have all Cartesian type and no other cell type
    AssertDimension(n_cell_types[0], n_cell_batches);

    // should have as many different Jacobians as we have cell batches in x
    // direction; as the mesh is a cube, we can easily calculate it
    deallog << "Number of different Jacobians: "
            << mf.get_mapping_info().cell_data[0].jacobians[0].size() / 2
            << std::endl;
  }

  // check again, now using a mapping that displaces points
  {
    MappingQ<dim> mapping(3);
    mf.reinit(mapping, dof, constraints, quad, data);

    std::vector<unsigned int> n_cell_types(4, 0);
    for (unsigned int i = 0; i < n_cell_batches; ++i)
      n_cell_types[mf.get_mapping_info().get_cell_type(i)]++;

    // should have all general type and no other cell type
    AssertDimension(n_cell_types[3], n_cell_batches);

    // should have as many different Jacobians as we have cell batches in x
    // direction times the number of quadrature points; as the mesh is a cube,
    // we can easily calculate it
    deallog << "Number of different Jacobians: "
            << mf.get_mapping_info().cell_data[0].jacobians[0].size()
            << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  deallog << std::setprecision(3);

  {
    deallog.push("2d");
    test<2>();
    test_cube<2>();
    test_parallelogram<2>();
    test_deformed_cube<2>();
    deallog.pop();
    deallog.push("3d");
    test<3>();
    test_cube<3>();
    test_parallelogram<3>();
    test_deformed_cube<3>();
    deallog.pop();
  }
}
