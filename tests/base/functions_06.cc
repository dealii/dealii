// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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


// Test VectorFunctionFromScalarFunctionObject

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>


template <int dim>
void check1 ()
{
  VectorFunctionFromScalarFunctionObject<dim>
  object (&Point<dim>::norm, 1, 3);

  Assert (object.n_components == 3, ExcInternalError());

  for (unsigned int i=0; i<10; ++i)
    {
      Point<dim> p;
      for (unsigned int d=0; d<dim; ++d)
        p[d] = i+d;

      for (unsigned int c=0; c<3; ++c)
        Assert (object.value(p,c) == (c==1 ? p.norm() : 0),
                ExcInternalError());

      Vector<double> v(3);
      object.vector_value (p, v);
      for (unsigned int c=0; c<3; ++c)
        Assert (v(c) == (c==1 ? p.norm() : 0),
                ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


template <int dim>
void check2 ()
{
  Point<dim> q;
  for (unsigned int d=0; d<dim; ++d)
    q[d] = d;

  ScalarFunctionFromFunctionObject<dim>
  object (std_cxx11::bind (&Point<dim>::distance,
                           q,
                           std_cxx11::_1));

  for (unsigned int i=0; i<10; ++i)
    {
      Point<dim> p;
      for (unsigned int d=0; d<dim; ++d)
        p[d] = i+d;

      Assert (object.value(p) == q.distance (p),
              ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int main()
{
  std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check1<1> ();
  check1<2> ();
  check1<3> ();

  check2<1> ();
  check2<2> ();
  check2<3> ();
}


