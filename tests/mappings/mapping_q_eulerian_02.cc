// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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


// test that DataOut and MappingQEulerian agree on how to output a
// displaced mesh. This was broken between r20158 and r21072
//
// this is a variant of _01 except that we use a quadratic FE field

#include <deal.II/base/function.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"



template <int dim>
class Displacement : public Function<dim>
{
public:
  Displacement()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return p[component] * p[0] / 2;
  }

  void
  vector_value(const Point<dim> &p, Vector<double> &v) const
  {
    for (unsigned int i = 0; i < dim; ++i)
      v(i) = p[i] * p[0] / 2;
  }
};


template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);

  FESystem<dim>   fe(FE_Q<dim>(2), dim);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> displacements(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler, Displacement<dim>(), displacements);

  MappingQEulerian<dim> euler(2, dof_handler, displacements);
  // now the actual test
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(displacements, "displacement");

  // output with all cells curved
  data_out.build_patches(euler, 1, DataOut<dim>::curved_inner_cells);
  data_out.write_gnuplot(deallog.get_file_stream());
}



int
main()
{
  initlog();
  deallog << std::setprecision(4);
  deallog.get_file_stream() << std::setprecision(4);

  test<1>();
  test<2>();
  test<3>();
}
