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


// test that DataOut and MappingQEulerian agree on how to output a
// displaced mesh. This was broken between r20158 and r21072

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/identity_matrix.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/multithread_info.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <vector>

using namespace dealii;


template <int dim>
class Displacement : public Function<dim>
{
public:
  Displacement() :
    Function<dim>(dim)
  {}

  double value (const Point<dim> &p,
                const unsigned int component) const
  {
    return p[component];
  }

  void vector_value (const Point<dim> &p,
                     Vector<double> &v) const
  {
    for (unsigned int i=0; i<dim; ++i)
      v(i) = p[i];
  }
};


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);

  FESystem<dim> fe(FE_Q<dim>(1),dim);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> displacements (dof_handler.n_dofs());

  VectorTools::interpolate (dof_handler,
                            Displacement<dim>(),
                            displacements);

  MappingQEulerian<dim> euler(2, displacements, dof_handler);
  // now the actual test
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  std::vector<std::string> names (dim, "displacement");
  data_out.add_data_vector(displacements,names);

  // output with all cells curved
  data_out.build_patches(euler,1,DataOut<dim>::curved_inner_cells);
  data_out.write_gnuplot(deallog.get_file_stream());
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision (4);
  logfile << std::setprecision (4);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
}


