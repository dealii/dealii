// shapes.cc,v 1.18 2003/04/09 15:49:55 wolf Exp
// (c) Guido Kanschat
//
// Show the shape functions implemented.

#include "../tests.h"
#include "shapes.h"
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_FaceQ_shape_functions()
{
  MappingQ1<dim> m;
  
  FE_FaceQ<dim> q0(0);
  plot_face_shape_functions(m, q0, "FaceQ0", update_values);
  
  FE_FaceQ<dim> q1(1);
  plot_face_shape_functions(m, q1, "FaceQ1", update_values);
  
  FE_FaceQ<dim> q2(2);
  plot_face_shape_functions(m, q2, "FaceQ2", update_values);
  
  FE_FaceQ<dim> q3(3);
  plot_face_shape_functions(m, q3, "FaceQ3", update_values);
  
  FE_FaceQ<dim> q4(4);
  plot_face_shape_functions(m, q4, "FaceQ4", update_values);
}


int
main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  plot_FE_FaceQ_shape_functions<2>();
  plot_FE_FaceQ_shape_functions<3>();
  
  return 0;
}



