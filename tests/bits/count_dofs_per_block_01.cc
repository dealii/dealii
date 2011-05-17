//----------------------------  count_dofs_per_block_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  count_dofs_per_block_01.cc  ---------------------------

// the DoFTools::count_dofs_per_{component,block} functions resized the output
// array to fe.n_components or fe.n_blocks even if a grouping argument was
// given. this would appear wrong and can lead to all sorts of interesting
// behavior if not caught by an assertion early enough


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <string>


std::string output_file_name = "count_dofs_per_block_01/output";



void print (const std::vector<unsigned int> &v)
{  
  deallog << v.size();
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << v[i];
  deallog << std::endl;
}

  

template <int dim>
void
check ()
{
                                   // create tria and dofhandler
                                   // objects. set different boundary
                                   // and sub-domain ids
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global (1);
  for (int i=0; i<2; ++i)
    {
      tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement ();
    }

				   // Taylor-Hood element
  FESystem<dim> fe (FE_Q<dim>(1), dim, FE_DGQ<dim>(0), 1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

				   // no grouping
  {
    std::vector<unsigned int> dpc(dim+1);
    DoFTools::count_dofs_per_component (dof_handler, dpc);
    print (dpc);
  }

  {
    std::vector<unsigned int> dpc(dim+1);
    DoFTools::count_dofs_per_block (dof_handler, dpc);
    print (dpc);
  }


				   // grouping into less groups than
				   // components
  {
    std::vector<unsigned int> group(dim+1, 0);
    group[dim] = 1;
    std::vector<unsigned int> dpc(2);
    DoFTools::count_dofs_per_component (dof_handler, dpc, false, group);
    Assert (dpc.size() == 2, ExcInternalError());
    print (dpc);
  }
  
  {
    std::vector<unsigned int> group(dim+1, 0);
    group[dim] = 1;
    std::vector<unsigned int> dpc(2);
    DoFTools::count_dofs_per_block (dof_handler, dpc, group);
    Assert (dpc.size() == 2, ExcInternalError());
    print (dpc);
  }  

				   // grouping into more groups than
				   // components
  {
    std::vector<unsigned int> group(dim+1, 2*dim);
    group[dim] = 0;
    std::vector<unsigned int> dpc(2*dim+1);
    DoFTools::count_dofs_per_component (dof_handler, dpc, false, group);
    Assert (dpc.size() == 2*dim+1, ExcInternalError());
    print (dpc);
  }
  
  {
    std::vector<unsigned int> group(dim+1, 2*dim);
    group[dim] = 0;
    std::vector<unsigned int> dpc(2*dim+1);
    DoFTools::count_dofs_per_block (dof_handler, dpc, group);
    Assert (dpc.size() == 2*dim+1, ExcInternalError());
    print (dpc);
  }  
}



int main()
{
  std::ofstream logfile(output_file_name.c_str());
  logfile << std::setprecision (2);
  deallog << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}
