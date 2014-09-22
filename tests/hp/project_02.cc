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


// check VectorTools::project for BlockVector<double> arguments and hp::


#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_accessor.h>

#include <fstream>


// define the multi-linear function x or x*y or x*y*z that we will
// subsequently project onto the ansatz space
template <int dim>
class F : public Function<dim>
{
public:
  virtual double value (const Point<dim> &p,
                        const unsigned int = 0) const
  {
    double s = 1;
    for (unsigned int i=0; i<dim; ++i)
      s *= p[i];
    return s;
  }
};


template<int dim>
void test()
{
  // create 2 triangulations with the
  // same coarse grid, and refine
  // them differently
  Triangulation<dim> tria;

  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));
  hp::DoFHandler<dim> dh (tria);
  dh.distribute_dofs (fe);

  BlockVector<double> v(2);
  v.block(0).reinit (dh.n_dofs()/2);
  v.block(1).reinit (dh.n_dofs()-dh.n_dofs()/2);
  v.collect_sizes();

  ConstraintMatrix cm;
  cm.close ();
  VectorTools::project (dh, cm, hp::QCollection<dim>(QGauss<dim>(3)), F<dim>(),
                        v);

  for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dh.begin_active();
       cell != dh.end(); ++cell)
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      {
        // check that the error is
        // somewhat small. it won't
        // be zero since we project
        // and do not interpolate
        Assert (std::fabs (v(cell->vertex_dof_index(i,0,cell->active_fe_index())) -
                           F<dim>().value (cell->vertex(i)))
                < 1e-4,
                ExcInternalError());
        deallog << cell->vertex(i) << ' ' << v(cell->vertex_dof_index(i,0,cell->active_fe_index()))
                << std::endl;
      }
}


int main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1>();
  test<2>();
  test<3>();
}

