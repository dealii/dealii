// $Id$
// (c) Guido Kanschat
//
// Check performance of finite elements

#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_lib.dg.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <vector>
#include <fstream>
#include <string>

template <int dim>
void performance (Triangulation<dim>& tr,
		  const FiniteElement<dim>& fe,
		  const Quadrature<dim>& quadrature,
		  UpdateFlags flags)
{
  deallog << "Create dofs" << endl;
  DoFHandler<dim> dof (tr);
  dof.distribute_dofs (fe);

  deallog << "Create FEValues" << endl;
  
  FEValues<dim> val (fe, quadrature, flags);

  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  DoFHandler<dim>::active_cell_iterator end = dof.end();

  deallog << "Loop" << endl;
  
  for (;cell != end ; ++cell)
    val.reinit(cell);

  deallog << "End" << endl;
}

int main ()
{
  deallog.log_execution_time(true);
  Triangulation<2> tr;
  GridGenerator::hyper_ball (tr);
  tr.refine_global (9);
  
  QGauss<2> gauss (3);
  
  FEQ2<2> element;

  deallog.push("points");
  performance (tr, element, gauss, update_q_points);
  deallog.pop();
  
  deallog.push("values");
  performance (tr, element, gauss, update_values);
  deallog.pop();
  
  deallog.push("grads-");
  performance (tr, element, gauss, update_gradients);
  deallog.pop();
  
  deallog.push("2nd---");
  performance (tr, element, gauss, update_second_derivatives);
  deallog.pop();

  deallog.push("matrix");
  performance (tr, element, gauss, UpdateFlags (update_q_points
	       | update_JxW_values
	       | update_values
	       | update_gradients));
  deallog.pop();

  deallog.push("all---");
  performance (tr, element, gauss, UpdateFlags (update_q_points
	       | update_JxW_values
	       | update_values
	       | update_gradients
	       | update_second_derivatives));
  deallog.pop();
}
