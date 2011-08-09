// shapes.cc,v 1.18 2003/04/09 15:49:55 wolf Exp
// (c) Guido Kanschat
//
// Show the shape functions implemented.

#include "../tests.h"
#include "shapes.h"
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_FaceQ_shape_functions()
{
  MappingQ1<dim> m;
  
  FE_FaceQ<dim> q0(0);
  FE_FaceQ<dim> q1(1);
  FE_FaceQ<dim> q2(2);
  FE_FaceQ<dim> q3(3);
  FE_DGQ<dim> dg0(0);
  FE_DGQ<dim> dg1(1);
  FESystem<dim> sys0(q0,1);
  FESystem<dim> sys1(q1,1);
  FESystem<dim> sys00(q0,1,dg0,1);
  FESystem<dim> sys11(q1,1,dg1,1);
  plot_face_shape_functions(m, q0, "FaceQ0", update_values);
  plot_face_shape_functions(m, q1, "FaceQ1", update_values);
  plot_face_shape_functions(m, q2, "FaceQ2", update_values);
  plot_face_shape_functions(m, q3, "FaceQ3", update_values);
  plot_face_shape_functions(m, sys0, "System0", update_values);  
  plot_face_shape_functions(m, sys1, "System1", update_values);  
  plot_face_shape_functions(m, sys00, "System0-0", update_values);
  plot_face_shape_functions(m, sys11, "System1-1", update_values);
}


int
main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
//  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  plot_FE_FaceQ_shape_functions<2>();
//  plot_FE_FaceQ_shape_functions<3>();
  
  return 0;
}



