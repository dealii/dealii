// nedelec_3.cc,v 1.2 2003/04/09 15:49:55 wolf Exp
// (c) Wolfgang Bangerth
//
// Some more tests with the Nedelec element, this time checking
// whether embedding matrices work correctly

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_tools.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 2




template<int dim>
inline void
check ()
{
  deallog << dim << "-D" << std::endl;
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

                                   // work on a once-refined grid
  tr.refine_global (1);

  FE_Nedelec<dim> fe_ned(1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe_ned);

                                   // generate a function on the
                                   // coarse grid (with shape function
                                   // values equal to the index of the
                                   // line)
  Vector<double> coarse_grid_values (fe_ned.dofs_per_cell);
  for (unsigned int i=0; i<coarse_grid_values.size(); ++i)
    coarse_grid_values(i) = i;

                                   // then transfer this function
                                   // from the coarse grid to the
                                   // child cells.
  Vector<double> values (dof.n_dofs());
  for (typename DoFHandler<dim>::active_cell_iterator
         c = dof.begin_active();
       c!=dof.end(); ++c)
    {
      Vector<double> tmp (fe_ned.dofs_per_cell);
      fe_ned.get_prolongation_matrix (c->index()).vmult (tmp, coarse_grid_values);
      c->set_dof_values (tmp, values);
    };
  

                                   // then output these values at the
                                   // quadrature points of all cells
                                   // of the finer grid
  QTrapez<dim> quadrature;
  std::vector<Vector<double> > shape_values (quadrature.n_quadrature_points,
                                             Vector<double>(dim));
  FEValues<dim> fe(fe_ned, quadrature,
                   update_values|update_q_points);

  for (typename DoFHandler<dim>::active_cell_iterator
         c = dof.begin_active();
       c!=dof.end(); ++c)
    {
      deallog << "  CELL " << c << std::endl;
      fe.reinit(c);
      fe.get_function_values (values, shape_values);

      for (unsigned int q=0; q<quadrature.n_quadrature_points; ++q)
        {
          deallog << "    q=" << q
                  << ", xq=" << fe.quadrature_point(q)
                  << ", f=[";
          for (unsigned int d=0; d<dim; ++d)
            deallog << (d==0 ? "" : " ")
                    << shape_values[q](d);
          
          deallog << "]" << std::endl;
        };
                
      deallog << std::endl;
    }
}


int
main()
{
  std::ofstream logfile ("nedelec_3.output");
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  check<2>();
  check<3>();
  
  return 0;
}



