// $Id$
// Copyright (C) 2001 Michael Stadler, Wolfgang Bangerth
//

#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1_eulerian.h>
#include <fe/fe_values.h>
#include <fe/fe_lib.lagrange.h>
#include <vector>
#include <fstream>
#include <string>



template<int dim>
inline void
show_values(FiniteElement<dim>& fe,
	    const char* name)
{
  deallog.push (name);
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 2., 5.);

				   // shift one point of the cell
				   // somehow
  if (dim > 1)
    tr.begin_active()->vertex(2)(dim-1) += 1./sqrt(2.);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);


				   // construct the MappingQ1Eulerian
				   // object
  
  FESystem<dim> mapping_fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> flowfield_dof_handler(tr);
  flowfield_dof_handler.distribute_dofs(mapping_fe);
  Vector<double> map_points(flowfield_dof_handler.n_dofs());
  MappingQ1Eulerian<dim> mapping(map_points, flowfield_dof_handler);


  QGauss2<dim> quadrature_formula;
  
  FEValues<dim> fe_values(mapping, fe, quadrature_formula,
			  UpdateFlags(update_values |
				      update_JxW_values |
				      update_gradients |
				      update_second_derivatives));
  DoFHandler<dim>::cell_iterator c = dof.begin();
  fe_values.reinit(c);

  for (unsigned int k=0; k<quadrature_formula.n_quadrature_points; ++k)
    {
      deallog << quadrature_formula.point(k) << std::endl;
      deallog << "\tJxW " << fe_values.JxW(k);
	  
      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
	{
	  deallog << "\tValues " << fe_values.shape_value(i,k);
	  deallog << "\tGrad " << fe_values.shape_grad(i,k);
	  deallog << "\t2nd " << fe_values.shape_2nd_derivative(i,k);
	  deallog << std::endl;
	}
    }
  deallog.pop ();
}



template<int dim>
void show_values()
{
  FE_Q<dim> q1(1);
  show_values(q1, "DGQ1");

  FE_Q<dim> q2(2);
  show_values(q2, "DGQ2");
}


int
main()
{
  std::ofstream logfile ("mapping_q1_eulerian.output");
  logfile.precision (2);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);



  deallog.push ("1d");
  show_values<1>();
  deallog.pop ();
  
  deallog.push ("2d");
  show_values<2>();
  deallog.pop ();
  
  deallog.push ("3d");
  show_values<3>();
  deallog.pop ();

  return 0;
}



