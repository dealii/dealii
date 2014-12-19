// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2014 by the deal.II authors
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



// like _27, but even simpler. this is based on a testcase proposed by
// Christoph Heiniger

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>


class F : public Function<2>
{
public:
  F() : Function<2>(2) {}

  virtual void vector_value (const Point<2> &p,
			     Vector<double> &v) const
    {
      // make the function equal to (0,x^2)
      v[0] = 0;
      v[1] = p[0]*p[0];
    }
};



Tensor<1,1> curl (const Tensor<2,2> &grads)
{
  return Point<1>(grads[1][0] - grads[0][1]);
}


Tensor<1,3> curl (const Tensor<2,3> &grads)
{
  return Point<3>(grads[2][1] - grads[1][2],
		  grads[0][2] - grads[2][0],
		  grads[1][0] - grads[0][1]);
}



template<int dim>
void test (const Triangulation<dim> &tr,
           const FiniteElement<dim> &fe)
{
  deallog << "FE=" << fe.get_name()
          << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  Vector<double> fe_function(dof.n_dofs());
  // set the elements of the vector in such a way that the function
  // equals the vector function (0,x^2)
  ConstraintMatrix cm;
  cm.close ();
  VectorTools::project (dof, cm, QGauss<2>(4), F(), fe_function);

  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature,
                           update_values | update_gradients | update_q_points);
  fe_values.reinit (dof.begin_active());

  // let the FEValues object compute the
  // divergences at quadrature points
  std::vector<typename dealii::internal::CurlType<dim>::type> curls (quadrature.size());
  std::vector<Tensor<2,dim> > grads (quadrature.size());
  FEValuesExtractors::Vector extractor(0);
  fe_values[extractor].get_function_curls (fe_function, curls);
  fe_values[extractor].get_function_gradients (fe_function, grads);

  // now compare
  for (unsigned int q=0; q<quadrature.size(); ++q)
    {
      deallog << "  curls[q]= " << curls[q] << std::endl
	      << "  grads[q]= " << grads[q] << std::endl;
      Assert ((curl(grads[q]) - curls[q]).norm()
	      <= 1e-10,
	      ExcInternalError());

      // we know the function F=(0,x^2) and we chose an element with
      // high enough degree to exactly represent this function, so the
      // curl of this function should be
      //   curl F = d_x F_y - d_y F_x = 2x
      // verify that this is true
      Assert (std::fabs(curls[q][0] - 2*fe_values.quadrature_point(q)[0])
	      <= 1e-10,
	      ExcInternalError());
    }
}



template<int dim>
void test_hyper_cube()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);

  FESystem<dim> fe(FE_Q<dim>(2), dim);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test_hyper_cube<2>();
}
