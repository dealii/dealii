// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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



/* Compare the kelly estimator for the same square but in codimension 0 and the other in codimension 1,
   Both should return the same estimator */



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>

template <int dim, int spacedim>
void
check ()
{
  Functions::CosineFunction<spacedim> function;
  Triangulation<dim,spacedim> tr;

  GridGenerator::hyper_cube(tr, -1,1);

  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();


  FE_Q<dim,spacedim> element(3);
  DoFHandler<dim,spacedim> dof(tr);
  dof.distribute_dofs(element);

  MappingQ<dim,spacedim> mapping(3);
  QGauss<dim-1> q_face(4);

  std::map<types::boundary_id,const Function<spacedim>*> neumann_bc;
  neumann_bc[0] = &function;

  Vector<double> v (dof.n_dofs());
  VectorTools::interpolate (mapping, dof, function, v);

  Vector<float> error (tr.n_active_cells());

  KellyErrorEstimator<dim,spacedim>::estimate (mapping, dof, q_face, neumann_bc,
                                               v, error);

  deallog << "Estimated error:" << std::endl;
  for (unsigned int i=0; i<error.size(); ++i)
    deallog << error(i)*100 << std::endl;
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);


  deallog.push ("2d_2");
  check<2,2> ();
  deallog.pop ();
  deallog.push ("2d_3");
  check<2,3> ();
  deallog.pop ();
}
