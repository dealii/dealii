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
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/mapping_cartesian.h>
#include <fe/mapping_q1.h>
#include <fe/mapping_q.h>
#include <fe/fe_values.h>
#include <vector>
#include <fstream>
#include <strstream>
#include <string>

template <int dim>
void performance (Triangulation<dim>& tr,
		  const Mapping<dim>& mapping,
		  const FiniteElement<dim>& fe,
		  const Quadrature<dim>& quadrature,
		  UpdateFlags flags)
{
  deallog << "Create dofs" << std::endl;
  DoFHandler<dim> dof (tr);
  dof.distribute_dofs (fe);

  deallog << "Create FEValues" << std::endl;
  
  FEValues<dim> val (mapping, fe, quadrature, flags);

  DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
  DoFHandler<dim>::active_cell_iterator end = dof.end();

  deallog << "Loop" << std::endl;
  
  for (;cell != end ; ++cell)
    val.reinit(cell);

  deallog << "End" << std::endl;
}

template <int dim>
void loop (std::vector<FiniteElement<dim> *> elements,
	    const Mapping<dim>& mapping,
	    Triangulation<dim>& tr)
{
  QGauss<dim> gauss (4);
  
  typename std::vector<FiniteElement<dim> *>::iterator elementp = elements.begin ();
  typename std::vector<FiniteElement<dim> *>::iterator end = elements.end ();

  for (;elementp != end; ++elementp)
    {
      const FiniteElement<dim>& element = **elementp;

      char dofs[20];
      std::ostrstream ost (dofs, 19);
      ost << element.n_dofs_per_cell() << std::ends;
      
      deallog.push(dofs);
      
      deallog.push("points");
      performance (tr, mapping, element, gauss, update_q_points);
      deallog.pop();
  
      deallog.push("values");
      performance (tr, mapping, element, gauss, update_values);
      deallog.pop();
  
      deallog.push("grads-");
      performance (tr, mapping, element, gauss, update_gradients);
      deallog.pop();
  
      deallog.push("2nd---");
      performance (tr, mapping, element, gauss, update_second_derivatives);
      deallog.pop();
  
      deallog.push("matrix");
      performance (tr, mapping, element, gauss, update_q_points
		   | update_JxW_values
		   | update_values
		   | update_gradients);
      deallog.pop();

      deallog.push("all---");
      performance (tr, mapping, element, gauss, update_q_points
		   | update_JxW_values
		   | update_values
		   | update_gradients
		   | update_second_derivatives);
      deallog.pop();
      deallog.pop();
    }
}

int main ()
{
  std::ofstream of ("performance.log");
  deallog.attach (of);
  deallog.log_execution_time(true);
  deallog.log_time_differences(true);
  Triangulation<2> tr;
  GridGenerator::hyper_ball (tr);
  tr.refine_global (8);
  
  MappingCartesian<2> cartesian;
  MappingQ1<2> mapping;
  MappingQ<2> mappingq1(1);
  MappingQ<2> mappingq2(2);
  std::vector<FiniteElement<2>*> el2d;
  el2d.push_back (new FE_Q<2> (1));
  el2d.push_back (new FE_Q<2> (2));
  el2d.push_back (new FE_Q<2> (3));
  el2d.push_back (new FE_Q<2> (4));

  loop (el2d, cartesian, tr);
  loop (el2d, mapping, tr);
  loop (el2d, mappingq1, tr);
  loop (el2d, mappingq2, tr);
}
