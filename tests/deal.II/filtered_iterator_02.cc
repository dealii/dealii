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
bool level_equal_to_3 (const Iterator c)
{
  return (static_cast<unsigned int>(c->level()) == 3);
}



template <typename Iterator>
bool level_equal_to (const Iterator     c,
                     const unsigned int level)
{
  return (static_cast<unsigned int>(c->level()) == level);
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
      const IteratorFilters::LevelEqualTo predicate(3);
      FilteredIterator<active_cell_iterator>
      begin = make_filtered_iterator(tria.begin_active(), predicate),
      end = make_filtered_iterator (static_cast<active_cell_iterator>(tria.end()), predicate);

      Assert (std::distance (begin, end) ==
              static_cast<signed int>(tria.n_active_cells (3)),
              ExcInternalError());
      logfile << "Check 1: "
              << (std::distance (begin, end) ==
                  static_cast<signed int>(tria.n_active_cells (3))
                  ?
                  "OK" : "Failed")
              << std::endl;
    };


  // check 2: count number of cells
  // on some level in a different way
  if (true)
    {
      bool (*predicate) (const active_cell_iterator)
        = &level_equal_to_3<active_cell_iterator>;
      FilteredIterator<active_cell_iterator>
      begin (predicate, tria.begin_active (3)),
            end   (predicate, tria.end());

      Assert (std::distance (begin, end) ==
              static_cast<signed int>(tria.n_active_cells (3)),
              ExcInternalError());
      logfile << "Check 2: "
              << (std::distance (begin, end) ==
                  static_cast<signed int>(tria.n_active_cells (3))
                  ?
                  "OK" : "Failed")
              << std::endl;
    };


  // check 3: count number of cells
  // on some level in yet a different
  // way
  if (true)
    {
      bool (*predicate) (const active_cell_iterator, const unsigned int)
        = &level_equal_to<active_cell_iterator>;
      FilteredIterator<active_cell_iterator>
      begin (std::bind2nd (std::ptr_fun(predicate), 3),
             tria.begin_active (3)),
                               end   (std::bind2nd (std::ptr_fun(predicate), 3),
                                      tria.end());

      Assert (std::distance (begin, end) ==
              static_cast<signed int>(tria.n_active_cells (3)),
              ExcInternalError());
      logfile << "Check 3: "
              << (std::distance (begin, end) ==
                  static_cast<signed int>(tria.n_active_cells (3))
                  ?
                  "OK" : "Failed")
              << std::endl;
    };


  // check 4: and yet another possibility
  if (true)
    {
      typedef FilteredIterator<active_cell_iterator> FI;

      bool (*predicate) (const active_cell_iterator, const unsigned int)
        = &level_equal_to<active_cell_iterator>;
      Assert (std::distance (FI(std::bind2nd (std::ptr_fun(predicate), 3))
                             .set_to_next_positive(tria.begin_active()),
                             FI(std::bind2nd (std::ptr_fun(predicate), 3), tria.end())) ==
              static_cast<signed int>(tria.n_active_cells (3)),
              ExcInternalError());
      logfile << "Check 4: "
              << (std::distance (FI(std::bind2nd (std::ptr_fun(predicate), 3))
                                 .set_to_next_positive(tria.begin_active()),
                                 FI(std::bind2nd (std::ptr_fun(predicate), 3), tria.end())) ==
                  static_cast<signed int>(tria.n_active_cells (3))
                  ?
                  "OK" : "Failed")
              << std::endl;
    };


  // check 5: check that we loop over
  // all cells with a given subdomain
  // id
  if (true)
    {
      typedef FilteredIterator<active_cell_iterator> FI;
      const IteratorFilters::SubdomainEqualTo predicate(1);
      FI cell (predicate);
      cell.set_to_next_positive (tria.begin_active());
      active_cell_iterator endc (tria.end());
      active_cell_iterator cell1 = tria.begin_active ();

      while (cell1->subdomain_id () != 1)
        ++cell1;

      while (true)
        {
          // move filtered iterator ahead
          ++cell;
          // move unfiltered iterator
          // ahead
          ++cell1;
          while ((cell1 != endc) &&
                 (cell1->subdomain_id () != 1))
            ++cell1;

          Assert (cell == cell1, ExcInternalError());
          Assert (cell1 == cell, ExcInternalError());

          if (cell.state() != IteratorState::valid)
            break;
        };
      Assert (cell == endc, ExcInternalError());
      Assert (cell1 == endc, ExcInternalError());

      logfile << "Check 5: OK" << std::endl;
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

