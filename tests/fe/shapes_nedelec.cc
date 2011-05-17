// shapes.cc,v 1.18 2003/04/09 15:49:55 wolf Exp
// (c) Guido Kanschat
//
// Show the shape functions implemented.

#include "../tests.h"
#include "shapes.h"
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_Nedelec_shape_functions()
{
  MappingQ1<dim> m;
  FE_Nedelec<dim> p0(0);
//   plot_shape_functions(m, p1, "Nedelec1");
//   plot_face_shape_functions(m, p1, "Nedelec1");
  test_compute_functions(m, p0, "Nedelec0");
  FE_Nedelec<dim> p1(1);
  test_compute_functions(m, p1, "Nedelec1");
}


int
main()
{
  std::ofstream logfile ("shapes_nedelec/output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  deallog << "FE_Nedelec<2>" << std::endl;
  plot_FE_Nedelec_shape_functions<2>();
  deallog << "FE_Nedelec<3>" << std::endl;
  plot_FE_Nedelec_shape_functions<3>();
  
  return 0;
}



