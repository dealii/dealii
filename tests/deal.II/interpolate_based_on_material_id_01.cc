// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// check VectorTools::interpolate_based_on_material_id

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_dgq.h>

#include <fstream>
#include <vector>


template <int dim>
class F :  public Function<dim>
{
public:
  F (const unsigned int q) : q(q) {}

  virtual double value (const Point<dim> &p,
                        const unsigned int) const
  {
    double v=0;
    for (unsigned int d=0; d<dim; ++d)
      for (unsigned int i=0; i<=q; ++i)
        v += (d+1)*(i+1)*std::pow (p[d], 1.*i);
    return v;
  }

private:
  const unsigned int q;
};



template <int dim>
void test ()
{
  Triangulation<dim>     triangulation;
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (3);
  std::map<types::material_id, const Function<dim>*> functions;
  for (typename Triangulation<dim>::active_cell_iterator
	 cell = triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    {
      cell->set_material_id(cell->index() % 128);
      if (functions.find(cell->index() % 128) == functions.end())
	functions[cell->index() % 128]
	  = new ConstantFunction<dim>(cell->index() % 128);
    }
  
  for (unsigned int p=1; p<7-dim; ++p)
    {
      FE_DGQ<dim>              fe(p);
      DoFHandler<dim>        dof_handler(triangulation);
      dof_handler.distribute_dofs (fe);

      Vector<double> interpolant (dof_handler.n_dofs());
      VectorTools::interpolate_based_on_material_id (MappingQ1<dim>(),
						     dof_handler,
						     functions,
						     interpolant);
      for (typename DoFHandler<dim>::active_cell_iterator
	     cell = dof_handler.begin_active();
	   cell != dof_handler.end();
	   ++cell)
	{
	  Vector<double> values (fe.dofs_per_cell);
	  cell->get_dof_values (interpolant, values);
	  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	    Assert (values[i] == cell->index() % 128, ExcInternalError());
	}
    }
  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
}

