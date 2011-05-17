//----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

// check MGTools::count_dofs_per_component

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_tools.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

void log_vector (const std::vector<std::vector<unsigned int> >& count)
{
  for (unsigned int l=0;l<count.size();++l)
    {
      deallog << "Level " << l;
      for (unsigned int c=0;c<count[l].size();++c)
	deallog << '\t' << count[l][c];
      deallog << std::endl;
    }
}

template <int dim>
void check_fe(FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement ();
  tr.refine_global(1);
  
  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  
  std::vector<std::vector<unsigned int> > count(tr.n_levels());
  MGTools::count_dofs_per_component(mgdof, count, false);
  log_vector(count);
  MGTools::count_dofs_per_component(mgdof, count, true);
  log_vector(count);

  std::vector<unsigned int> target(fe.n_components());
  for (unsigned int i=0;i<target.size();++i)
    target[i] = i/3;
  deallog << std::endl << "Target";
  for (unsigned int i=0;i<target.size();++i)
    deallog << '\t' << target[i];
  deallog << std::endl;
  
  MGTools::count_dofs_per_component(mgdof, count, false, target);
  log_vector(count);
  MGTools::count_dofs_per_component(mgdof, count, true, target);
  log_vector(count);
  
}


template <int dim>
void check()
{
  FE_Q<dim> q1(1);
  FE_Q<dim> q2(2);
  FE_DGQ<dim> dq1(1);
  
  FESystem<dim> s1(q1, 2, q2,1);

  check_fe(s1);
  if (dim>1)
    {
      FE_RaviartThomas<dim> rt(1);
      FESystem<dim> s2(rt, 2, dq1,1);
      FESystem<dim> s3(rt,1, s1, 2);
      
      check_fe(s2);
      check_fe(s3);
    }
}

int main()
{
  std::ofstream logfile("count_01/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}
