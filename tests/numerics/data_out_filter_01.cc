// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

// This test shows the ability of DataOutFilter to merge vertex duplicates taking
// into account floating point inprecision.
// The test creates a distorted triangulation with vertices very close to each other
// that are subsequently merged by the DataOutFilter class.
// Note that while the test case will likely never happen, there are cases
// in which vertices during mesh output have slightly different positions
// that would not be merged if the points would simply be compared for equality.

#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/numerics/data_out.h>




template <int dim>
void
test ()
{
  Triangulation<dim> tria(Triangulation<dim>::MeshSmoothing::none,true);
  GridGenerator::hyper_cube(tria, 1., 1. + 1000.0 * std::numeric_limits<double>::epsilon());

  FE_Q<dim> fe1(1);

  DoFHandler<dim> dof1(tria);
  dof1.distribute_dofs(fe1);

  Vector<double> v1(dof1.n_dofs());
  for (unsigned int i=0; i<v1.size(); ++i) v1(i) = i;

  DataOut<dim> data_out;
  data_out.add_data_vector (dof1, v1, "linear");
  data_out.build_patches ();

  DataOutBase::DataOutFilter data_filter
  (DataOutBase::DataOutFilterFlags (true, false));

  data_out.write_filtered_data (data_filter);

  deallog << "Number of filtered nodes: " << data_filter.n_nodes() << std::endl;

  Triangulation<dim> tria2(Triangulation<dim>::MeshSmoothing::none,true);
  GridGenerator::hyper_cube(tria2, 1., 1 + 10.0 * std::numeric_limits<double>::epsilon());

  FE_Q<dim> fe2(1);

  DoFHandler<dim> dof2(tria2);
  dof2.distribute_dofs(fe2);

  Vector<double> v2(dof2.n_dofs());
  for (unsigned int i=0; i<v2.size(); ++i) v2(i) = i;

  DataOut<dim> data_out2;
  data_out2.add_data_vector (dof2, v2, "linear");
  data_out2.build_patches ();

  DataOutBase::DataOutFilter data_filter2
  (DataOutBase::DataOutFilterFlags (true, false));

  data_out2.write_filtered_data (data_filter2);

  deallog << "Number of filtered nodes: " << data_filter2.n_nodes() << std::endl;

  deallog << "ok" << std::endl;
}


int main(int argc, char *argv[])
{
  initlog();
  test<2>();
  test<3>();

  return 0;
}
