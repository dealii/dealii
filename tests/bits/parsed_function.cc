// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// This program tests the functionality of the function parser
// wrapper.

#include "../tests.h"
#include <fstream>
#include <iomanip>
#include <map>
#include <deal.II/base/logstream.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

template <int dim>
void Test()
{
  // A parameter handler
  ParameterHandler prm;

  // A triangulation
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria,0,1);
  tria.refine_global(3);

  // Vertices
  const std::vector<Point<dim> > &vertices = tria.get_vertices();

  // Test vector declaration
  for (unsigned int i=0; i<dim; ++i)
    {
      std::string id =
        "Function " + Utilities::int_to_string(dim) +
        " - " + Utilities::int_to_string(i);
      prm.enter_subsection(id);

      Functions::ParsedFunction<dim>::declare_parameters(prm, i+1);
      prm.set("Function constants", "f=" + Utilities::int_to_string(i+1));

      // It is cos(pi f x_i) *t^i, numbering from zero
      std::string expr = "cos(pi*f*x)*t";
      if (i>0) expr+=  "; cos(pi*f*y)*t";
      if (i>1) expr+=  "; cos(pi*f*z)*t";

      prm.set("Function expression", expr);

      Functions::ParsedFunction<dim> function(i+1);
      function.parse_parameters(prm);

      prm.leave_subsection();

      // Now test the difference from t=0 to t=1
      for (double t=0.; t<1; t+= .1)
        {
          function.set_time(t);
          Point<dim> p;
          std::vector<Vector<double> > values(vertices.size(),
                                              Vector<double>(i+1));
          function.vector_value_list(vertices, values);
          for (unsigned int j=0; j<vertices.size(); ++j)
            {
              double delta=0.;
              for (unsigned int di=0; di<i; ++di)
                {
                  delta = values[j](di) - std::cos(numbers::PI*(i+1)*vertices[j][di])*t;
                  if (std::abs(delta) > 1e-10)
                    deallog << "p(" << di << "): " << vertices[j] << ", delta: "
                            << delta << std::endl;
                }
            }
        }
    }
  deallog << "Tested on: " << std::endl;
  prm.log_parameters(deallog);
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Test<1>();
  Test<2>();
  Test<3>();
}




