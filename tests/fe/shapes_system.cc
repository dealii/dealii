// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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


#include "../tests.h"
#include "shapes.h"
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2

template<int dim>
void plot_FE_System_shape_functions()
{
  MappingQ1<dim> m;

//   FESystem<dim> p1(FE_Q<dim>(2), 1,
//                    FE_Q<dim>(dim<3 ? 3 : 2), 2);
//   plot_shape_functions(m, p1, "System1");
//   plot_face_shape_functions(m, p1, "System1");
//   test_compute_functions(m, p1, "System1");

//   FESystem<dim> p2(FE_Q<dim>(2), 1,
//                    FESystem<dim> (FE_Q<dim>(1),1,
//                                   FE_DGP<dim>(3),3,
//                                   FE_DGQ<dim>(0),2), 2,
//                    FE_DGQ<dim>(0), 2);
//   plot_shape_functions(m, p2, "System2");
//   plot_face_shape_functions(m, p2, "System2");
//   test_compute_functions(m, p2, "System2");

  // some tests with the Nedelec
  // element. don't try to make sense
  // out of the composed elements,
  // they are simply constructed as
  // complicated as possible to
  // trigger as many assertions as
  // possible (and they _have_, in
  // the past, literally dozens of
  // assertions)
  if (dim != 1)
    {
      FESystem<dim> p3(FE_Nedelec<dim>(1), 1,
                       FESystem<dim> (FE_Q<dim>(1),1,
                                      FE_DGP<dim>(3),3,
                                      FE_Nedelec<dim>(1),2), 2,
                       FE_DGQ<dim>(0), 2);
      test_compute_functions(m, p3, "System_Nedelec_1");

      // the following is simply too
      // expensive in 3d...
      if (dim != 3)
        {
          FESystem<dim> p4(p3, 1,
                           FESystem<dim> (FE_Q<dim>(1),1,
                                          p3,3,
                                          FE_Nedelec<dim>(1),2), 1,
                           p3, 1);
          test_compute_functions(m, p4, "System_Nedelec_2");
        };
    };
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  deallog << "FE_System<1>" << std::endl;
  plot_FE_System_shape_functions<1>();
  deallog << "FE_System<2>" << std::endl;
  plot_FE_System_shape_functions<2>();
  deallog << "FE_System<3>" << std::endl;
  plot_FE_System_shape_functions<3>();

  return 0;
}



