// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2015 by the deal.II authors
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



// call VectorTools::create_boundary_rhs. this crashed at the time of writing
// this because we forgot to resize the size of a local_dofs array on each
// cell. test and fix by Daniel Gerecht, University of Heidelberg


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <fstream>
#include <numeric>


void test ()
{
  //make grid
  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  //define DoFhandler and FEs
  FE_Q<2> u(2),u2(QIterated<1>(QTrapez<1>(),3));

  hp::FECollection<2> fe_collection;
  fe_collection.push_back(u);
  fe_collection.push_back(u2);

  hp::FECollection<2> fe_collection2;
  fe_collection2.push_back(u);

  hp::DoFHandler<2> hp_dof_handler(triangulation);
  hp::DoFHandler<2> hp_dof_handler2(triangulation);

  //set different fe for testing
  hp::DoFHandler<2>::active_cell_iterator
  cell=hp_dof_handler.begin_active(),endc = hp_dof_handler.end();

  for (; cell!=endc; ++cell)
    cell->set_active_fe_index(1);
  hp_dof_handler.begin_active()->set_active_fe_index(0);

  //distribute dofs
  hp_dof_handler.distribute_dofs(fe_collection);
  hp_dof_handler2.distribute_dofs(fe_collection2);

  // define input for create_rhs function
  Vector<double> rhs_vector(hp_dof_handler.n_dofs());
  Vector<double> rhs_vector2(hp_dof_handler2.n_dofs());

  types::boundary_id myints[1] = {0};
  std::set<types::boundary_id> boundary_ids  (myints,myints+1);
  myints[0]=0;
  hp::QCollection<1>   quadrature;
  quadrature.push_back (QGauss<1>(1));

  // choose rhs=1. because the circumference
  // of the domain is four and the shape
  // functions all add up to one, if you take
  // the sum over the elements of the
  // resulting rhs vector, you need to get
  // four
  ConstantFunction< 2 > rhs_function(1);

  VectorTools::create_boundary_right_hand_side ( hp_dof_handler2,
                                                 quadrature,
                                                 rhs_function,
                                                 rhs_vector2,
                                                 boundary_ids);
  Assert (std::fabs (std::accumulate (rhs_vector2.begin(),
                                      rhs_vector2.end(), 0.)
                     - 4) < 1e-12,
          ExcInternalError());
  deallog << rhs_vector2.l2_norm() << std::endl;

  VectorTools::create_boundary_right_hand_side ( hp_dof_handler,
                                                 quadrature,
                                                 rhs_function,
                                                 rhs_vector,
                                                 boundary_ids);
  Assert (std::fabs (std::accumulate (rhs_vector.begin(),
                                      rhs_vector.end(), 0.)
                     - 4) < 1e-12,
          ExcInternalError());
  deallog << rhs_vector.l2_norm() << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  logfile.precision(2);
  deallog << std::setprecision(2);

  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
