// shapes.cc,v 1.18 2003/04/09 15:49:55 wolf Exp
// (c) Guido Kanschat
//
// Show the shape functions implemented.

#include "../tests.h"
#include "shapes.h"
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_Q_shape_functions()
{
  MappingQ1<dim> m;

  FE_Q_iso_Q1<dim> q1(1);
  plot_shape_functions(m, q1, "Q1");
  plot_face_shape_functions(m, q1, "Q1");
  test_compute_functions(m, q1, "Q1");

  FE_Q_iso_Q1<dim> q2(2);
  plot_shape_functions(m, q2, "Q2_iso_Q1");
  plot_face_shape_functions(m, q2, "Q2_iso_Q1");
  test_compute_functions(m, q2, "Q2_iso_Q1");

				   // skip the following tests to
				   // reduce run-time
  if (dim < 3)
    {
      FE_Q_iso_Q1<dim> q3(3);
      plot_shape_functions(m, q3, "Q3_iso_Q1");
      plot_face_shape_functions(m, q3, "Q3_iso_Q1");
      test_compute_functions(m, q3, "Q3_iso_Q1");

      FE_Q_iso_Q1<dim> q4(4);
      plot_shape_functions(m, q4, "Q4_iso_Q1");
      plot_face_shape_functions(m, q4, "Q4_iso_Q1");
      test_compute_functions(m, q4, "Q4_iso_Q1");
    };

}


int
main()
{
  std::ofstream logfile ("shapes_q_iso_q1/output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  plot_FE_Q_shape_functions<1>();
  plot_FE_Q_shape_functions<2>();
  plot_FE_Q_shape_functions<3>();

  return 0;
}



