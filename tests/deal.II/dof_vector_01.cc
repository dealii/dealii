// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

// Test constructor and access of DoFVector as well as sync for Vector
// and BlockVector

#include "../tests.h"
#include <deal.II/base/mg_level_object.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/dofs/dof_vector.h>

#include <cmath>
#include <cstdlib>
#include <fstream>

template <class DH, typename number>
void test(const DH& dh)
{
  deallog << dh.get_fe().get_name() << std::endl;
  
  DoFVector<DH, Vector<number> > dv1(dh);
  dv1.sync();
  Vector<number>& v1 = dv1.data();
  deallog << "dim " << DH::dimension << " v1: " << v1.size() << std::endl;
  
  const DoFVector<DH, Vector<number> > dv2(dh);
  const Vector<number>& v2 = dv2.data();
  
  Vector<number> w(5);
  const DoFVector<DH, Vector<number> > dv3(dh, w);
  const Vector<number>& v3 = dv3.data();

  DoFVector<DH, BlockVector<number> > dvb(dh);
  dvb.sync();
  BlockVector<number>& vb = dvb.data();
  deallog << "dim " << DH::dimension << " vb: " << vb.get_block_indices().to_string() << std::endl;
  
}

template <class DH, typename number>
void test_mg(const DH& dh)
{
  deallog << "MG" << dh.get_fe().get_name() << std::endl;
  
  DoFVector<DH, MGLevelObject<Vector<number> > > dv1(dh);
  dv1.sync();
  MGLevelObject<Vector<number> >& v1 = dv1.data();
  deallog << "dim " << DH::dimension << " v1: (" << v1.min_level()
	  << "->" << v1.max_level() << ')' << std::endl;
  for (unsigned int i=v1.min_level(); i <= v1.max_level();++i)
    deallog << " l" << i << ": " << v1[i].size() << std::endl;
  
  const DoFVector<DH, Vector<number> > dv2(dh);
  const Vector<number>& v2 = dv2.data();
  
  Vector<number> w(5);
  const DoFVector<DH, Vector<number> > dv3(dh, w);
  const Vector<number>& v3 = dv3.data();

  DoFVector<DH, MGLevelObject<BlockVector<number> > > dvb(dh);
  dvb.sync();
  MGLevelObject<BlockVector<number> >& vb = dvb.data();
  deallog << "dim " << DH::dimension << " vb: (" << vb.min_level()
	  << "->" << vb.max_level() << ')' << std::endl;
  for (unsigned int i=vb.min_level(); i <= vb.max_level();++i)
    deallog << " l" << i << ": " << vb[i].get_block_indices().to_string() << std::endl;
}

template <int dim>
void test_grid()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  FE_Q<dim> fe(1);
  
  DoFHandler<dim> dh(tr);
  dh.distribute_dofs(fe);
  dh.distribute_mg_dofs(fe);
  
  test<DoFHandler<dim>, float>(dh);
  test_mg<DoFHandler<dim>, float>(dh);
  test<DoFHandler<dim>, double>(dh);
  test_mg<DoFHandler<dim>, double>(dh);

  FESystem<dim> fes(fe,3);
  DoFHandler<dim> ds(tr);
  ds.distribute_dofs(fes);
  ds.distribute_mg_dofs(fes);
  
  test<DoFHandler<dim>, float>(ds);
  test_mg<DoFHandler<dim>, float>(ds);
  test<DoFHandler<dim>, double>(ds);
  test_mg<DoFHandler<dim>, double>(ds);
}


int main ()
{
  initlog();
  test_grid<2>();
  test_grid<3>();
  return 0;
}

