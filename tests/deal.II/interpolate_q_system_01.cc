// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
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



// check that VectorTools::interpolate works for FE_System(FE_Q(p)) elements correctly on
// a uniformly refined mesh for functions of degree q

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
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <vector>


template <int dim>
class F :  public Function<dim>
{
public:
  F (const unsigned int q) : Function<dim>(3), q(q) {}

  virtual void vector_value (const Point<dim> &p,
                             Vector<double>   &v) const
  {
    for (unsigned int c=0; c<v.size(); ++c)
      {
        v(c) = 0;
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int i=0; i<=q; ++i)
            v(c) += (d+1)*(i+1)*std::pow (p[d], 1.*i)+c;
      }
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

  for (unsigned int p=1; p<6-dim; ++p)
    {
      FE_Q<dim> fe_1(p);
      FE_Q<dim> fe_2(p+1);
      FESystem<dim>              fe(fe_1, 2,
                                    fe_2, 1);
      DoFHandler<dim>        dof_handler(triangulation);
      dof_handler.distribute_dofs (fe);

      Vector<double> interpolant (dof_handler.n_dofs());
      Vector<float>  error (triangulation.n_active_cells());
      for (unsigned int q=0; q<=p+2; ++q)
        {
          // interpolate the function
          VectorTools::interpolate (dof_handler,
                                    F<dim> (q),
                                    interpolant);

          // then compute the interpolation error
          VectorTools::integrate_difference (dof_handler,
                                             interpolant,
                                             F<dim> (q),
                                             error,
                                             QGauss<dim>(q+2),
                                             VectorTools::L2_norm);
          if (q<=p)
            Assert (error.l2_norm() < 1e-12*interpolant.l2_norm(),
                    ExcInternalError());

          deallog << fe.get_name() << ", P_" << q
                  << ", rel. error=" << error.l2_norm() / interpolant.l2_norm()
                  << std::endl;
        }
    }
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

