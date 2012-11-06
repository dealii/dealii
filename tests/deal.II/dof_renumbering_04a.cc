//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2007, 2008, 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// an extract of the _04 test that produced different (integer)
// results on different machines, inexplicably


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/multigrid/mg_dof_accessor.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>



template <int dim>
void
print_dofs (const DoFHandler<dim> &dof)
{
  std::vector<unsigned int> v (dof.get_fe().dofs_per_cell);

  for (typename DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
       cell != dof.end(); ++cell)
    {
      cell->get_dof_indices (v);
      deallog << "cell=" << cell << std::endl;
      for (unsigned int i=0; i<v.size(); ++i)
	deallog << v[i] << std::endl;
    }
}




template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, -1., 1.);
  tr.refine_global (1);

  FESystem<dim> fe(FE_DGQ<dim>(1),1);
  DoFHandler<dim> dof(tr);
  
  dof.distribute_dofs(fe);
  DoFRenumbering::boost::Cuthill_McKee(dof);
  print_dofs (dof);
}


int main ()
{
  std::ofstream logfile ("dof_renumbering_04a/output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  check<2> ();
}
