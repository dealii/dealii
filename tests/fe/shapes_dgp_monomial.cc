// shapes.cc,v 1.18 2003/04/09 15:49:55 wolf Exp
// (c) Guido Kanschat
//
// Show the shape functions implemented.

#include "../tests.h"
#include "shapes.h"
#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_DGPMonomial_shape_functions()
{
  MappingQ1<dim> m;

  FE_DGPMonomial<dim> p1(1);
  plot_shape_functions(m, p1, "DGPMonomial1");
  plot_face_shape_functions(m, p1, "DGPMonomial1");
  test_compute_functions(m, p1, "DGPMonomial1");

  FE_DGPMonomial<dim> p2(2);
  plot_shape_functions(m, p2, "DGPMonomial2");
  plot_face_shape_functions(m, p2, "DGPMonomial2");
  test_compute_functions(m, p2, "DGPMonomial2");

  if (dim<3)
    {
      FE_DGPMonomial<dim> p3(3);
      plot_shape_functions(m, p3, "DGPMonomial3");
      plot_face_shape_functions(m, p3, "DGPMonomial3");
      test_compute_functions(m, p3, "DGPMonomial3");
    }
}


int
main()
{
  std::ofstream logfile ("shapes_dgp_monomial/output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  plot_FE_DGPMonomial_shape_functions<1>();
  plot_FE_DGPMonomial_shape_functions<2>();
  plot_FE_DGPMonomial_shape_functions<3>();
  
  return 0;
}



