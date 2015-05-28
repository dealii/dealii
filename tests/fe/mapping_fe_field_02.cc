// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2015 by the deal.II authors
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

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>

#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#define PRECISION 2

// We create two FESystem with spacedim components with an FE_Bernstein
// and an FE_Q and test if the two position vectors are the same for a
// simple mesh like an hyper_cube.

template<int dim, int spacedim>
void test(const unsigned int degree) {
  deallog << "dim = " << dim << ", spacedim = " << spacedim << std::endl;
  deallog << "degree = " << degree << std::endl;

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global (2);

  FE_Bernstein<dim,spacedim> fe(degree);
  FE_Q<dim,spacedim> feq(degree);
  FESystem<dim,spacedim> fe_sys(fe, spacedim);
  FESystem<dim,spacedim> fe_sysq(feq, spacedim);

  DoFHandler<dim,spacedim> dof_sys(tria);
  DoFHandler<dim,spacedim> dof_sysq(tria);

  dof_sys.distribute_dofs(fe_sys);
  dof_sysq.distribute_dofs(fe_sysq);

  Vector<double> euler;
  Vector<double> eulerq;

  euler.reinit(dof_sys.n_dofs());
  eulerq.reinit(dof_sysq.n_dofs());

  const ComponentMask mask(spacedim, true);

  VectorTools::get_position_vector(dof_sys, euler, mask);
  VectorTools::get_position_vector(dof_sysq, eulerq, mask);
  
  for (unsigned int i=0; i<euler.size(); ++i)
    Assert (std::fabs(euler[i] - eulerq[i]) < 1e-12, ExcInternalError());

  deallog << "OK" << std::endl;
}

int main()
{
  initlog();

  for(unsigned int d=1; d<5; ++d) {
    test<1,1>(d);
    test<1,2>(d);
    test<2,2>(d);
    test<2,3>(d);
    test<3,3>(d);
  }  
}
