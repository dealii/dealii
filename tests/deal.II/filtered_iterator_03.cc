// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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


// check filtered iterators


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/filtered_iterator.h>

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>


std::ofstream logfile("output");


DeclException2 (ExcNumberMismatch,
                int, int,
                << "The numbers " << arg1 << " and " << arg2
                << " should be equation, but are not.");



typedef Triangulation<2>::active_cell_iterator active_cell_iterator;

template <typename Iterator>
bool always_true (const Iterator)
{
  return true;
}





void test ()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global (1);
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global (2);

  // we now have a number of cells,
  // flag them with some subdomain
  // ids based on their position, in
  // particular we take the quadrant
  // (octant)
  active_cell_iterator cell = tria.begin_active (),
                       endc = tria.end ();
  for (; cell!=endc; ++cell)
    {
      unsigned int subdomain = 0;
      for (unsigned int d=0; d<2; ++d)
        if (cell->center()(d) > 0)
          subdomain |= (1<<d);
      Assert (subdomain < (1<<2), ExcInternalError());

      cell->set_subdomain_id (subdomain);
    };


  // check 1: count number of cells
  // on some level
  if (true)
    {
      FilteredIterator<active_cell_iterator>
      begin = make_filtered_iterator(tria.begin_active(), &always_true<active_cell_iterator>),
      end = make_filtered_iterator (static_cast<active_cell_iterator>(tria.end()), &always_true<active_cell_iterator>);

      Assert (std::distance (begin, end) ==
              static_cast<signed int>(tria.n_active_cells ()),
              ExcInternalError());
      deallog << std::distance(begin,end) << ' '
              << tria.n_active_cells() << std::endl;
      logfile << "Check 1: "
              << (std::distance (begin, end) ==
                  static_cast<signed int>(tria.n_active_cells ())
                  ?
                  "OK" : "Failed")
              << std::endl;
    };
}


int main ()
{
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  return 0;
}

