// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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



// check VectorTools::create_right_hand_side, but this time for a scalar
// primitive element


#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/numerics/vector_tools.h>



template <int dim>
class MySquareFunction : public Function<dim>
{
public:
  MySquareFunction () : Function<dim>(1) {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const
  {
    return (component+1)*p.square();
  }

  virtual void   vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const
  {
    values(0) = value(p,0);
  }
};




template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  if (dim==2)
    {
      GridGenerator::hyper_ball(tr, Point<dim>(), 1);
      tr.reset_manifold(0);
    }
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  // create a system element composed
  // of one Q1 and one Q2 element
  FE_Q<dim> element(1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  // use a more complicated mapping
  // of the domain and a quadrature
  // formula suited to the elements
  // we have here
  MappingQ<dim> mapping(3);

  QGauss<dim> quadrature(3);

  Vector<double> rhs (dof.n_dofs());
  VectorTools::create_right_hand_side (dof, quadrature,
                                       MySquareFunction<dim>(),
                                       rhs);
  for (unsigned int i=0; i<rhs.size(); ++i)
    deallog << rhs(i) << std::endl;
}



int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (8);
  deallog << std::fixed;
  deallog.attach(logfile);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
