//----------------------------------------------------------------------------
//    transfer.cc,v 1.13 2005/12/30 16:07:03 guido Exp
//    Version:
//
//    Copyright (C) 2000 - 2007, 2009, 2010, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

//TODO:[GK] Add checks for RT again!

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

template <int dim>
void check_simple(const FiniteElement<dim>& fe, std::ostream& os)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  tr.begin(2)->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);

  MGTransferPrebuilt<Vector<double> > transfer;
  transfer.build_matrices(mgdof);
  
  transfer.print_matrices(os);
}


int main()
{
  std::ofstream logfile("transfer_01/output");
  deallog << std::setprecision(6);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check_simple (FE_DGP<2>(0), logfile);
  check_simple (FE_DGP<2>(1), logfile);
  check_simple (FE_DGQ<2>(1), logfile);
  check_simple (FE_DGQ<2>(2), logfile);
  check_simple (FE_Q<2>(1), logfile);
  check_simple (FE_Q<2>(2), logfile);
  check_simple (FESystem<2>(FE_DGQ<2>(1), 2), logfile);
  check_simple (FESystem<2>(FE_DGP<2>(1), 2, FE_DGQ<2>(1), 3), logfile);

  check_simple (FE_RaviartThomasNodal<2>(1), logfile);
//  check_simple (FESystem<2>(FE_RaviartThomas<2>(1),1,FE_DGQ<2>(0),2), logfile);

  check_simple (FE_DGQ<3>(1), logfile);
  check_simple (FE_Q<3>(2), logfile);
}
