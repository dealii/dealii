//----------------------------  rtn_2.cc  ---------------------------
//    rtn_2.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2003, 2005, 2006, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rtn_2.cc  ---------------------------

// rt(2) had some problems with shape functions (tested in rt_7.cc). because
// it is so simple, just run the same test for the RT-Nodal element as well

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 2



template<int dim>
void
plot_shape_functions(const unsigned int degree)
{
  FE_RaviartThomasNodal<dim> fe_rt(degree);
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

  DoFHandler<dim> dof(tr);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(fe_rt);
      
  QTrapez<1> q_trapez;
  const unsigned int div=10;
  QIterated<dim> q(q_trapez, div);
  FEValues<dim> fe(fe_rt, q, update_values|update_gradients|update_q_points);
  fe.reinit(c);

  for (unsigned int q_point=0; q_point<q.size(); ++q_point)
    {
      if (q_point % QIterated<1>(q_trapez,div).size() == 0)
        deallog << std::endl;
      
      deallog << fe.quadrature_point(q_point) << " ";
	      
      for (unsigned int i=0;i<fe_rt.dofs_per_cell;++i)
        for (unsigned int c=0; c<fe.get_fe().n_components(); ++c)
          deallog << " " << fe.shape_value_component(i,q_point,c);

      deallog << std::endl;
    }
}


int
main()
{
  std::ofstream logfile ("rtn_2/output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  plot_shape_functions<2>(2);
  
  return 0;
}



