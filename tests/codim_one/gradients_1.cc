// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2015 by the deal.II authors
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



// Controls that the covariant matrix is calculated properly. It uses
// a Q1 finite element to calculate the scalar product of the gradient
// of a projected function (a monomial) with the tangential to the
// cell surface taken in the cell midpoint.  The result obtained is
// compared with the exact one in the <1,2> case.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <string>

// all include files needed for the program

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <cmath>



std::ofstream logfile("output");

template <int dim, int spacedim>
void test(std::string filename, unsigned int degree = 1)

{
  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim> gi;

  gi.attach_triangulation (triangulation);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  // finite elements used for the
  // projection
  const FE_Q<dim,spacedim> fe (QIterated<1>(QTrapez<1>(),degree));
  const MappingQ<dim, spacedim> mapping(degree);

  DoFHandler<dim,spacedim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);

  deallog
      << "no. of cells "<< triangulation.n_cells() <<std::endl;
  deallog
      << "no. of dofs "<< dof_handler.n_dofs()<< std::endl;
  deallog
      << "no. of dofs per cell "<< fe.dofs_per_cell<< std::endl;


  //  definition of the exact function
  //  and calculation of the projected
  //  one
  Vector<double> projected_one(dof_handler.n_dofs());

  Functions::CosineFunction<spacedim> the_function;

  // Tensor<1,spacedim> exp;
  // exp[0]=1;
  // exp[1]=0;
  // if(spacedim==3)
  //     exp[2]=0;
  // Functions::Monomial<spacedim> the_function(exp);

  const QGauss<dim> quad(2*fe.degree+1);
  ConstraintMatrix constraints;
  constraints.close();
  VectorTools::project(mapping, dof_handler, constraints,
                       quad, the_function, projected_one);

  deallog << "L2 norm of projected vector: "
          << projected_one.l2_norm() << std::endl;


  // compute the H1 difference
  Vector<float> difference_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference (dof_handler, projected_one,
                                     the_function, difference_per_cell,
                                     quad, VectorTools::H1_norm);

  deallog << "H1 error: " << difference_per_cell.l2_norm() << std::endl;
}



int main ()
{
  logfile.precision (4);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-12);

  deallog<<"Test <1,2>, Q1, Q2, Q3"<<std::endl;
  test<1,2>(SOURCE_DIR "/grids/circle_4.inp",1);
  test<1,2>(SOURCE_DIR "/grids/circle_4.inp",2);
  test<1,2>(SOURCE_DIR "/grids/circle_4.inp",3);

  deallog<<std::endl;

  deallog<<"Test <2,3>, Q1, Q2, Q3"<<std::endl;
  test<2,3>(SOURCE_DIR "/grids/sphere_1.inp",1);
  test<2,3>(SOURCE_DIR "/grids/sphere_1.inp",2);
  test<2,3>(SOURCE_DIR "/grids/sphere_1.inp",3);


  return 0;
}
