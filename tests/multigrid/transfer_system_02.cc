//----------------------------------------------------------------------------
//    transfer.cc,v 1.13 2005/12/30 16:07:03 guido Exp
//    Version:
//
//    Copyright (C) 2000 - 2007, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------------

// like _01, but with a TransferSelect that selects all blocks

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_component.h>
#include <deal.II/multigrid/mg_tools.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

template <int dim, typename number, int spacedim>
void
reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
	       MGLevelObject<dealii::Vector<number> > &v)
{
  for (unsigned int level=v.get_minlevel();
       level<=v.get_maxlevel();++level)
    {
      unsigned int n = mg_dof.n_dofs (level);
      v[level].reinit(n);
    }

}


template <typename Transfer>
void
make_matrix (const Transfer &transfer,
	     const unsigned int high_level,
	     FullMatrix<double> &matrix)
{
  Vector<double> src (matrix.n());
  Vector<double> dst (matrix.m());
  for (unsigned int i=0; i<src.size(); ++i)
    {
      src = 0;
      src(i) = 1;
      transfer.prolongate (high_level, dst, src);
      for (unsigned int j=0; j<dst.size(); ++j)
	matrix(j,i) = dst(j);
    }
}



void print_matrix (const FullMatrix<double> &m)
{
  for (unsigned int i=0; i<m.m(); ++i)
    {
      for (unsigned int j=0; j<m.n(); ++j)
	deallog << m(i,j) << ' ';
      deallog << std::endl;
    }
}




template <int dim>
void check (const FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  MGDoFHandler<dim> mg_dof_handler(tr);
  mg_dof_handler.distribute_dofs(fe);

  DoFRenumbering::component_wise (mg_dof_handler);
  for (unsigned int level=0; level<tr.n_levels(); ++level)
    DoFRenumbering::component_wise (mg_dof_handler, level);

  MGTransferSelect<double> transfer;

				   // group all components into one
				   // block, and do the transfer
				   // matrices on this block
  transfer.build_matrices(mg_dof_handler, mg_dof_handler,
			  0, 0,
			  std::vector<unsigned int>(2, 0U),
			  std::vector<unsigned int>(2, 0U));

  FullMatrix<double> prolong_0_1 (mg_dof_handler.n_dofs(1),
				  mg_dof_handler.n_dofs(0));
  FullMatrix<double> prolong_1_2 (mg_dof_handler.n_dofs(2),
				  mg_dof_handler.n_dofs(1));

  deallog << "Level 0->1" << std::endl;
  make_matrix (transfer, 1, prolong_0_1);
  print_matrix (prolong_0_1);

  deallog << std::endl;

  deallog << "Level 1->2" << std::endl;
  make_matrix (transfer, 2, prolong_1_2);
  print_matrix (prolong_1_2);
}


int main()
{
  std::ofstream logfile("transfer_system_02/output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

//TODO: do in 1d
  check (FESystem<2>(FE_Q<2>(1),2));
}
