// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// check MGDoFAccessor::get_mg_dof_indices

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/function_map.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_tools.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;


template <int dim>
void dofs(const MGDoFHandler<dim> &dof)
{
  typename MGDoFHandler<dim>::cell_iterator cell;
  const typename MGDoFHandler<dim>::cell_iterator end = dof.end();

  std::vector<types::global_dof_index> indices;

  for (cell = dof.begin(); cell != end; ++cell)
    {
      indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_mg_dof_indices(indices);

      deallog << "Level " << cell->level();
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        deallog << " v" << cell->vertex(i);
      deallog << " dofs ";
      for (unsigned int i=0; i<indices.size(); ++i)
        deallog << ' ' << indices[i];
      deallog << std::endl;
    }
}


template <int dim>
void check_fe(FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  ZeroFunction<dim> zero;
  typename FunctionMap<dim>::type fmap;
  fmap.insert(std::make_pair(0, &zero));

  MGDoFHandler<dim> mgdof(tr);
  mgdof.distribute_dofs(fe);
  dofs(mgdof);
}


template <int dim>
void check()
{
  FE_Q<dim> q1(1);
  FE_Q<dim> q2(2);
//  FE_DGQ<dim> dq1(1);

  FESystem<dim> s1(q1, 2, q2,1);

  check_fe(q1);
  check_fe(q2);
  check_fe(s1);
}

int main()
{
  initlog(__FILE__);
  check<1> ();
  check<2> ();
  check<3> ();
}
