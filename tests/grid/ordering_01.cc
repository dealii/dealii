// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2013 by the deal.II authors
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

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <cmath>

using namespace dealii;


template <int dim>
void
show_ordering(const Triangulation<dim> &tr)
{
  for (typename Triangulation<dim>::cell_iterator cell = tr.begin(); cell != tr.end(); ++cell)
    for (typename Triangulation<dim>::cell_iterator other = tr.begin(); other != tr.end(); ++other)
      {
        deallog << (cell < other ? "less " : "not  ");
        deallog << cell->level_subdomain_id() << ':' << other->level_subdomain_id() << ' ';
        if (cell->active())
          deallog << cell->subdomain_id();
        else
          deallog << 'X';
        deallog << ':';
        if (other->active())
          deallog << other->subdomain_id();
        else
          deallog << 'X';
        deallog << ' ';
        deallog << cell->level() << ':' << other->level() << ' ';
        deallog << cell->index() << ':' << other->index() << ' ';
        deallog << std::endl;
      }
}

template <int dim>
void test1()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  typename Triangulation<dim>::active_cell_iterator cell = tr.begin_active();
  cell->set_subdomain_id(4);
  cell->set_level_subdomain_id(4);
  ++cell;
  cell->set_subdomain_id(3);
  cell->set_level_subdomain_id(5);

  tr.refine_global(1);
  cell = tr.begin_active();
  cell->set_level_subdomain_id(3);
  (++cell)->set_level_subdomain_id(4);
  (++cell)->set_level_subdomain_id(5);
  ++cell;
  (++cell)->set_level_subdomain_id(3);
  (++cell)->set_level_subdomain_id(4);
  (++cell)->set_level_subdomain_id(5);

  show_ordering(tr);
}

int main()
{
  initlog();

  test1<2>();
  return 0;
}
