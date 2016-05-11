// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2015 by the deal.II authors
//
// Author: Jihuan Tian <jihuan_tian@hotmail.com>
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


// Test curl related functions in integrators/maxwell.h

#include <deal.II/base/logstream.h>
#include <cstring>
#include <fstream>
#include <iostream>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/integrators/maxwell.h>

using namespace dealii;
using namespace LocalIntegrators::Maxwell;

template <int dim>
void make_grid(Triangulation<dim> &tr)
{
  GridGenerator::subdivided_hyper_cube(tr, 1, -1., 1.);
}

template <int dim>
void TestMaxwellCurl(Triangulation<dim> &tr)
{
  int element_order = 1;
  int quadrature_order = 3;

  DoFHandler<dim> dof_handler(tr);
  FESystem<dim> fe(FE_Q<dim>(element_order), dim);
  FEValues<dim> fe_values(fe, QGauss<dim>(quadrature_order), update_values | update_JxW_values | update_gradients | update_quadrature_points);
  FEFaceValues<dim> fe_face_values(fe, QGauss<dim-1>(quadrature_order), update_values | update_JxW_values | update_gradients | update_normal_vectors | update_quadrature_points);

  dof_handler.distribute_dofs(fe);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  FullMatrix<double> curl_curl_check(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> curl_check(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> nitsche_curl_check(dofs_per_cell, dofs_per_cell);

  curl_curl_check = 0;
  curl_check = 0;
  nitsche_curl_check = 0;

  typename::DoFHandler<dim>::cell_iterator cell = dof_handler.begin(0);

  fe_values.reinit(cell);

  curl_curl_matrix<dim>(curl_curl_check, fe_values, 1.);
  deallog << "curl_curl_matrix" << std::endl;
  curl_curl_check.print(deallog, 10);

  curl_matrix(curl_check, fe_values, fe_values, 1.);
  deallog << "curl_matrix" << std::endl;
  curl_check.print(deallog, 10);

  for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
    {
      fe_face_values.reinit(cell, face);
      nitsche_curl_matrix<dim>(nitsche_curl_check, fe_face_values, face, 1., 1.);
    }

  deallog << "nitsche_curl_matrix" << std::endl;
  nitsche_curl_check.print(deallog, 10);

  dof_handler.clear();
}

int main()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  Triangulation<3> tr;
  make_grid(tr);
  TestMaxwellCurl(tr);
}
