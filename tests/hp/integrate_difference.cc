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



// call VectorTools::integrate_difference with fe's distributed in the
// same random way as in hp/random


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>



template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (2);
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global (4-dim);

  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim> q_collection;
  for (unsigned int i=1; i<=4; ++i)
    {
      fe_collection.push_back(FE_Q<dim> (QIterated<1>(QTrapez<1>(),i)));
      q_collection.push_back (QGauss<dim> (i+2));
    }


  hp::DoFHandler<dim> dof_handler(tria);

  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    cell->set_active_fe_index (Testing::rand() % fe_collection.size());

  dof_handler.distribute_dofs(fe_collection);

  Vector<double> vec (dof_handler.n_dofs());
  for (unsigned int i=0; i<vec.size(); ++i)
    vec(i) = i;

  Vector<float> diff (tria.n_active_cells());

  VectorTools::NormType norms[] =
  {
    VectorTools::mean,
    VectorTools::L1_norm,
    VectorTools::L2_norm,
    VectorTools::Linfty_norm,
    VectorTools::H1_seminorm,
    VectorTools::W1p_seminorm
  };
  for (unsigned int i=0; i<sizeof(norms)/sizeof(norms[0]); ++i)
    {
      VectorTools::integrate_difference (dof_handler,
                                         vec,
                                         Functions::SquareFunction<dim>(),
                                         diff,
                                         q_collection,
                                         norms[i]);
      deallog << "i=" << i << ", diff=" << diff.l2_norm() << std::endl;
    }
}


int main ()
{
  std::ofstream logfile("output");
  logfile.precision(2);
  deallog << std::setprecision(2);

  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
