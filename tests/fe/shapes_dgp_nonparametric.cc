// shapes.cc,v 1.18 2003/04/09 15:49:55 wolf Exp
// (c) Guido Kanschat
//
// Show the shape functions implemented.

#include "../tests.h"
#include "shapes.h"
#include <fe/fe_dgp_nonparametric.h>
#include <fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_DGPNonparametric_shape_functions()
{
  MappingQ1<dim> m;

  FE_DGPNonparametric<dim> p0(0);
  plot_shape_functions(m, p0, "DGPNonparametric0");
  plot_face_shape_functions(m, p0, "DGPNonparametric0");

  FE_DGPNonparametric<dim> p1(1);
  plot_shape_functions(m, p1, "DGPNonparametric1");
  plot_face_shape_functions(m, p1, "DGPNonparametric1");

  FE_DGPNonparametric<dim> p2(2);
  plot_shape_functions(m, p2, "DGPNonparametric2");
  plot_face_shape_functions(m, p2, "DGPNonparametric2");
      
//    FE_DGPNonparametric<dim> p3(3);
//    plot_shape_functions(m, p3, "DGPNonparametric3");
//    plot_face_shape_functions(m, p3, "DGPNonparametric3");
      
//    FE_DGPNonparametric<dim> p4(4);
//    plot_shape_functions(m, p4, "DGPNonparametric4");
//    plot_face_shape_functions(m, p4, "DGPNonparametric4");
}


int
main()
{
  std::ofstream logfile ("shapes_dgp_nonparametric/output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  plot_FE_DGPNonparametric_shape_functions<1>();
  plot_FE_DGPNonparametric_shape_functions<2>();
  plot_FE_DGPNonparametric_shape_functions<3>();
  
  return 0;
}



