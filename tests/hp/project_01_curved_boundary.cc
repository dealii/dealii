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


// check VectorTools::project for Vector<double> arguments and hp


#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_accessor.h>

#include <fstream>


// define a function that is zero within the square
//   [-sqrt(2)/2,sqrt(2)/2]^d
// which is the domain produced by GridGenerator::hyper_sphere
// when using a linear mapping. we will want to see a solution
// that is not zero when using a higher order mapping and a Q2
// element
template <int dim>
class F : public Function<dim>
{
public:
  virtual double value (const Point<dim> &p,
                        const unsigned int = 0) const
  {
    // compute the linfty norm of p
    double m = 0;
    for (unsigned int d=0; d<dim; ++d)
      m = std::max(m, std::fabs(p[d]));

    // let the value increase linearly away from the square
    return std::max (0., m-std::sqrt(2.)/2);
  }
};


template<int dim>
void test()
{
  HyperBallBoundary<dim> boundary;

  // create 2 triangulations with the
  // same coarse grid, and refine
  // them differently
  Triangulation<dim> tria;
  GridGenerator::hyper_ball (tria);
  tria.set_boundary (0, boundary);


  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));
  hp::DoFHandler<dim> dh (tria);
  dh.distribute_dofs (fe);

  Vector<double> v(dh.n_dofs());

  ConstraintMatrix cm;
  cm.close ();

  // use the implicit Q1 mapping. this will yield a zero solution
  {
    VectorTools::project (dh, cm, hp::QCollection<dim>(QGauss<dim>(3)), F<dim>(),
                          v);
    deallog << v.l2_norm() << std::endl;
    Assert (v.l2_norm() == 0, ExcInternalError());
  }

  // use an explicit Q1 mapping. this will yield a zero solution
  {
    VectorTools::project (hp::MappingCollection<dim>(MappingQ1<dim>()),
                          dh, cm, hp::QCollection<dim>(QGauss<dim>(3)), F<dim>(),
                          v);
    deallog << v.l2_norm() << std::endl;
    Assert (v.l2_norm() == 0, ExcInternalError());
  }

  // use an explicit Q2 mapping. this will yield a nonzero solution if it's a
  // straight projection since some of the quadrature points are lying outside
  // the area where the function is zero
  {
    VectorTools::project (hp::MappingCollection<dim>(MappingQ<dim>(2)),
                          dh, cm, hp::QCollection<dim>(QGauss<dim>(3)), F<dim>(),
                          v);
    deallog << v.l2_norm() << std::endl;
    Assert (v.l2_norm() != 0, ExcInternalError());
  }

  // use an explicit Q2 mapping but enforce zero boundary values. this will
  // yield a nonzero solution since some of the quadrature points are lying
  // outside the area where the function is zero, even though the values of
  // the DoFs at the boundary are in fact zero (they are interpolated only at
  // points where the function is zero)
  {
    VectorTools::project (hp::MappingCollection<dim>(MappingQ<dim>(2)),
                          dh, cm, hp::QCollection<dim>(QGauss<dim>(3)), F<dim>(),
                          v, true);
    deallog << v.l2_norm() << std::endl;
    Assert (v.l2_norm() != 0, ExcInternalError());
  }

  // use an explicit Q2 mapping and project to the boundary first. this will
  // yield a nonzero solution since some of the quadrature points are lying
  // outside the area where the function is zero. furthermore, the values on
  // the boundary will be nonzero since we project, even though the
  // *interpolation* of boundary values onto the trace of the Q1 space is zero
  {
    VectorTools::project (hp::MappingCollection<dim>(MappingQ<dim>(2)),
                          dh, cm, hp::QCollection<dim>(QGauss<dim>(3)), F<dim>(),
                          v, false,
                          hp::QCollection<dim-1>(QGauss<dim-1>(2)), true);
    deallog << v.l2_norm() << std::endl;
    Assert (v.l2_norm() != 0, ExcInternalError());
    for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dh.begin_active();
         cell != dh.end(); ++cell)
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        deallog << cell->vertex(i) << ' ' << v(cell->vertex_dof_index(i,0,cell->active_fe_index()))
                << std::endl;
  }


  // same as above, but use a projection with a QTrapez formula. this happens
  // to evaluate the function only at points where it is zero, and
  // consequently the values at the boundary should be zero
  {
    VectorTools::project (hp::MappingCollection<dim>(MappingQ<dim>(2)),
                          dh, cm, hp::QCollection<dim>(QGauss<dim>(3)), F<dim>(),
                          v, false,
                          hp::QCollection<dim-1>(QTrapez<dim-1>()), true);
    deallog << v.l2_norm() << std::endl;
    Assert (v.l2_norm() != 0, ExcInternalError());
    for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dh.begin_active();
         cell != dh.end(); ++cell)
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
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

  test<2>();
  test<3>();
}

