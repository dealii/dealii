// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// We used to forget to instantiate a few static const member variables.
// Verify that this is now fixed.


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <class ACCESSOR>
LogStream &
operator<<(LogStream &log, const TriaIterator<ACCESSOR> &i)
{
  log << ACCESSOR::dimension << ' ' << ACCESSOR::structure_dimension << ' '
      << ACCESSOR::space_dimension << ' ';
  i.print(log);
  return log;
}


template <int dim>
void
test_in_dim(const DoFHandler<dim> &d1, const DoFHandler<dim> &d2)
{
  typename DoFHandler<dim>::active_cell_iterator a = d1.begin_active();
  typename DoFHandler<dim>::cell_iterator        l =
    d1.begin(d1.get_triangulation().n_levels() - 1);

  deallog << "a " << a << std::endl << "l " << l << std::endl;
}


template <int dim>
void
init_tria(Triangulation<dim> &tr)
{
  GridGenerator::hyper_cube(tr);
  tr.refine_global(4 - dim);
}


int
main()
{
  initlog();

  Triangulation<2> t2;
  init_tria(t2);

  FE_Q<2>       q21(1);
  DoFHandler<2> d21(t2);
  d21.distribute_dofs(q21);

  FE_Q<2>       q22(2);
  DoFHandler<2> d22(t2);
  d22.distribute_dofs(q22);

  test_in_dim(d21, d22);
}
