// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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


// the same as point_difference_01, but for the alternative point_difference
// algorithm

// check that VectorTools::point_difference returns the same before and after
// the change to a logarithmic algorithm



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
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <cmath>
#include <iomanip>


template<int dim>
class MySquareFunction : public Function<dim>
{
public:
  MySquareFunction () : Function<dim> () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const
  {
    return (component+1)*p.square()+1;
  }

  virtual void   vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const
  {
    values(0) = value(p,0);
  }
};


template<int dim>
class MyExpFunction : public Function<dim>
{
public:
  MyExpFunction () : Function<dim> () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const
  {
    return std::exp (p(0));
  }

  virtual void   vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const
  {
    values(0) = value(p,0);
  }
};



template <int dim>
void make_mesh (Triangulation<dim> &tria)
{

  GridGenerator::hyper_cube(tria, -1, 1);

  // refine the mesh in a random way so as to
  // generate as many cells with
  // hanging nodes as possible
  tria.refine_global (4-dim);
  const double steps[4] = { /*d=0*/ 0, 7, 3, 3 };
  for (unsigned int i=0; i<steps[dim]; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = tria.begin_active();
      for (unsigned int index=0; cell != tria.end(); ++cell, ++index)
        if (index % (3*dim) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement ();
    }
}




template <int dim>
void
check ()
{
  MappingQ1<dim> mapping;
  Triangulation<dim> tria;
  make_mesh (tria);

  FE_Q<dim> element(3);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(element);

  // test with two different functions: one
  // that is exactly representable on the
  // chosen finite element space, and one
  // that isn't
  for (unsigned int i=0; i<2; ++i)
    {
      static const MySquareFunction<dim>          function_1;
      static const Functions::CosineFunction<dim> function_2;

      const Function<dim> &
      function = (i==0 ?
                  static_cast<const Function<dim>&>(function_1) :
                  static_cast<const Function<dim>&>(function_2));

      Vector<double> v (dof.n_dofs());
      VectorTools::interpolate (dof, function, v);

      // for the following points, check the
      // function value, output it, and
      // verify that the value retrieved from
      // the interpolated function is close
      // enough to that of the real function
      //
      // also verify that the actual value is
      // roughly correct
      Point<dim> p[3];
      for (unsigned int d=0; d<dim; ++d)
        {
          p[0][d] = 0;
          p[1][d] = 0.5;
          p[2][d] = 1./3.;
        }
      Vector<double> difference(1);
      for (unsigned int i=0; i<3; ++i)
        {
          VectorTools::point_difference (mapping, dof, v, function, difference, p[i]);
          deallog << difference(0) << std::endl;
          Assert (difference(0) < 1e-4, ExcInternalError());

          VectorTools::point_difference (mapping, dof, v, ZeroFunction<dim>(),
                                         difference, p[i]);
          deallog << difference(0) << std::endl;
          Assert (std::abs(-difference(0) - function.value(p[i])) < 1e-4,
                  ExcInternalError());
        }
    }

  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (4);
  deallog.attach(logfile);
  deallog.depth_console (0);

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
