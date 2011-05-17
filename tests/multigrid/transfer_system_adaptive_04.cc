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

// test like _03 but with boundary conditions 
#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
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


template<int dim>
void refine_mesh (Triangulation<dim> &triangulation)
{
    bool cell_refined = false;
  for (typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active();
      cell != triangulation.end(); ++cell)
  {
      const Point<dim> p = cell->center();
      bool positive = p(0) > 0;
      if (positive)
      {
        cell->set_refine_flag();
        cell_refined = true;
      }
  }
  if(!cell_refined)//if no cell was selected for refinement, refine global
    for (typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active();
        cell != triangulation.end(); ++cell)
      cell->set_refine_flag();
  triangulation.execute_coarsening_and_refinement ();
}



template <int dim>
void check (const FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;

  std::vector<unsigned int> subdivisions (dim, 1);
  subdivisions[0] = 2;

  const Point<dim> bottom_left = (dim == 2 ?
      Point<dim>(-1,-1) : Point<dim>(-1,-1,-1));
  const Point<dim> top_right   = (dim == 2 ?
      Point<dim>(1,1) : Point<dim>(1,1,1));
  GridGenerator::subdivided_hyper_rectangle (tr,
      subdivisions, bottom_left, top_right, true);
  refine_mesh(tr);

  MGDoFHandler<dim> mg_dof_handler(tr);
  mg_dof_handler.distribute_dofs(fe);

  deallog << "Global  dofs: " << mg_dof_handler.n_dofs() << std::endl;
  for(unsigned int l=0; l<tr.n_levels(); ++l)
    {
      deallog << "Level " << l << " dofs:";
	deallog << ' ' << mg_dof_handler.n_dofs(l);
      deallog << std::endl;
    }

  DoFRenumbering::component_wise (mg_dof_handler);
  for (unsigned int level=0; level<tr.n_levels(); ++level)
    DoFRenumbering::component_wise (mg_dof_handler, level);

  std::vector<std::set<unsigned int> > boundary_indices(tr.n_levels());
  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    dirichlet_bc(fe.n_components());
  dirichlet_boundary[3] =             &dirichlet_bc;

  MGTools::make_boundary_list (mg_dof_handler, dirichlet_boundary,
			       boundary_indices);

  std::vector<unsigned int> block_selected(2,0U);
  MGTransferSelect<double> transfer;
  transfer.build_matrices(mg_dof_handler, 
      mg_dof_handler, 0,0, block_selected, 
      block_selected, boundary_indices);

  FullMatrix<double> prolong_0_1 (mg_dof_handler.n_dofs(1),
				  mg_dof_handler.n_dofs(0));

  deallog << "Level 0->1" << std::endl;
  make_matrix (transfer, 1, prolong_0_1);
  print_matrix (prolong_0_1);
}


int main()
{
  std::ofstream logfile("transfer_system_adaptive_04/output");
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

//TODO: do in 1d
  check (FESystem<2>(FE_Q<2>(1),2));
}
