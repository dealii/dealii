// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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

template <int dim, typename number, int spacedim>
void
reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
               MGLevelObject<dealii::Vector<number> > &v)
{
  for (unsigned int level=v.min_level();
       level<=v.max_level(); ++level)
    {
      unsigned int n = mg_dof.n_dofs (level);
      v[level].reinit(n);
    }

}


template <int dim>
void check_simple(const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);

  MGTransferPrebuilt<Vector<double> > transfer;
  transfer.build_matrices(mgdof);

  MGLevelObject<Vector<double> > u(0, tr.n_levels()-1);
  reinit_vector(mgdof, u);
  // First prolongate the constant
  // vector.  For Lagrange elements,
  // the values are just the number
  // of degrees of freedom.
  u[0] = 1;
  transfer.prolongate(1,u[1],u[0]);
  transfer.prolongate(2,u[2],u[1]);
  deallog << "u0\t" <<  (u[0]*u[0]+.5) << std::endl
          << "u1\t" <<  (u[1]*u[1]+.5) << std::endl
          << "u2\t" <<  (u[2]*u[2]+.5) << std::endl;
  // Now restrict the same vectors.
  u[1] = 0.;
  u[0] = 0.;
  transfer.restrict_and_add(2,u[1],u[2]);
  transfer.restrict_and_add(1,u[0],u[1]);
  deallog << "u1\t" <<  (u[1]*u[1]+.5) << std::endl
          << "u0\t" <<  (u[0]*u[0]+.5) << std::endl;

  // Now the same for a non-constant
  // vector
  for (unsigned int i=0; i<u[0].size(); ++i)
    u[0](i) = i;
  transfer.prolongate(1,u[1],u[0]);
  transfer.prolongate(2,u[2],u[1]);
  deallog << "u0\t" <<  (u[0]*u[0]+.5) << std::endl
          << "u1\t" <<  (u[1]*u[1]+.5) << std::endl
          << "u2\t" <<  (u[2]*u[2]+.5) << std::endl;
  // Now restrict the same vectors.
  u[1] = 0.;
  u[0] = 0.;
  transfer.restrict_and_add(2,u[1],u[2]);
  transfer.restrict_and_add(1,u[0],u[1]);
  deallog << "u1\t" <<  (u[1]*u[1]+.5) << std::endl
          << "u0\t" <<  (u[0]*u[0]+.5) << std::endl;

  // Fill a global vector by counting
  // from one up
  Vector<double> v;
  v.reinit (mgdof.n_dofs());
  for (unsigned int i=0; i<v.size(); ++i)
    v(i) = i+1;

  transfer.copy_to_mg(mgdof, u, v);
  for (unsigned int i=0; i<u[2].size(); ++i)
    deallog << ' ' << (int) u[2](i);
  deallog << std::endl;

  // Now do the opposite: fill a
  // multigrid vector counting the
  // dofs and see where the numbers go
  v = 0.;
  for (unsigned int i=0; i<u[2].size(); ++i)
    u[2](i) = i+1;
  transfer.copy_from_mg(mgdof, v, u);
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << ' ' << (int) v(i);
  deallog << std::endl;
  v.equ(-1., v);
  transfer.copy_from_mg_add(mgdof, v, u);
  deallog << "diff " << v.l2_norm() << std::endl;
}


int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(6);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check_simple (FE_DGP<2>(0));
  check_simple (FE_DGP<2>(1));
  check_simple (FE_DGQ<2>(1));
  check_simple (FE_DGQ<2>(2));
  check_simple (FE_Q<2>(1));
  check_simple (FE_Q<2>(2));
  check_simple (FESystem<2>(FE_DGQ<2>(1), 2));
  check_simple (FESystem<2>(FE_DGP<2>(1), 2, FE_DGQ<2>(1), 3));

  check_simple (FE_RaviartThomasNodal<2>(1));
//  check_simple (FESystem<2>(FE_RaviartThomas<2>(1),1,FE_DGQ<2>(0),2));

  check_simple (FE_DGQ<3>(1));
  check_simple (FE_Q<3>(2));
}
