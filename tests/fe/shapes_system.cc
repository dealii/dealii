// shapes.cc,v 1.18 2003/04/09 15:49:55 wolf Exp
// (c) Guido Kanschat
//
// Show the shape functions implemented.

#include "../tests.h"
#include "shapes.h"
#include <fe/fe_q.h>
#include <fe/fe_dgp.h>
#include <fe/fe_dgq.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1.h>
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
  std::ofstream logfile ("shapes_system/output");
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



