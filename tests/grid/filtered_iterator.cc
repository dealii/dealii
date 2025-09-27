// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check filtered iterators


#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <algorithm>
#include <numeric>

#include "../tests.h"


DeclException2(ExcNumberMismatch,
               int,
               int,
               << "The numbers " << arg1 << " and " << arg2
               << " should be equation, but are not.");



using active_cell_iterator = Triangulation<2>::active_cell_iterator;

template <typename Iterator>
bool
level_equal_to_3(const Iterator c)
{
  return (static_cast<unsigned int>(c->level()) == 3);
}



template <typename Iterator>
bool
level_equal_to(const Iterator c, const unsigned int level)
{
  return (static_cast<unsigned int>(c->level()) == level);
}


void
test()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(2);

  // we now have a number of cells,
  // flag them with some subdomain
  // ids based on their position, in
  // particular we take the quadrant
  // (octant)
  active_cell_iterator cell = tria.begin_active(), endc = tria.end();
  for (; cell != endc; ++cell)
    {
      unsigned int subdomain = 0;
      for (unsigned int d = 0; d < 2; ++d)
        if (cell->center()[d] > 0)
          subdomain |= (1 << d);
      AssertThrow(subdomain < (1 << 2), ExcInternalError());

      cell->set_subdomain_id(subdomain);
    };

  std::ostream &logfile = deallog.get_file_stream();

  // check 1: count number of cells
  // on some level
  if (true)
    {
      const IteratorFilters::LevelEqualTo    predicate(3);
      FilteredIterator<active_cell_iterator> begin(predicate),
        end(predicate, tria.end());
      begin.set_to_next_positive(tria.begin_active());

      Assert(std::distance(begin, end) ==
               static_cast<signed int>(tria.n_active_cells(3)),
             ExcInternalError());
      logfile << "Check 1: "
              << (std::distance(begin, end) ==
                      static_cast<signed int>(tria.n_active_cells(3)) ?
                    "OK" :
                    "Failed")
              << std::endl;
    };


  // check 2: count number of cells
  // on some level in a different way
  if (true)
    {
      bool (*predicate)(const active_cell_iterator) =
        &level_equal_to_3<active_cell_iterator>;
      FilteredIterator<active_cell_iterator> begin(predicate,
                                                   tria.begin_active(3)),
        end(predicate, tria.end());

      Assert(std::distance(begin, end) ==
               static_cast<signed int>(tria.n_active_cells(3)),
             ExcInternalError());
      logfile << "Check 2: "
              << (std::distance(begin, end) ==
                      static_cast<signed int>(tria.n_active_cells(3)) ?
                    "OK" :
                    "Failed")
              << std::endl;
    };


  // check 3: count number of cells
  // on some level in yet a different
  // way
  if (true)
    {
      bool (*predicate)(const active_cell_iterator, const unsigned int) =
        &level_equal_to<active_cell_iterator>;
      FilteredIterator<active_cell_iterator> begin(
        std::bind(predicate, std::placeholders::_1, 3), tria.begin_active(3)),
        end(std::bind(predicate, std::placeholders::_1, 3), tria.end());

      Assert(std::distance(begin, end) ==
               static_cast<signed int>(tria.n_active_cells(3)),
             ExcInternalError());
      logfile << "Check 3: "
              << (std::distance(begin, end) ==
                      static_cast<signed int>(tria.n_active_cells(3)) ?
                    "OK" :
                    "Failed")
              << std::endl;
    };


  // check 4: and yet another possibility
  if (true)
    {
      using FI = FilteredIterator<active_cell_iterator>;

      bool (*predicate)(const active_cell_iterator, const unsigned int) =
        &level_equal_to<active_cell_iterator>;
      Assert(std::distance(FI(std::bind(predicate, std::placeholders::_1, 3))
                             .set_to_next_positive(tria.begin_active()),
                           FI(std::bind(predicate, std::placeholders::_1, 3),
                              tria.end())) ==
               static_cast<signed int>(tria.n_active_cells(3)),
             ExcInternalError());
      logfile
        << "Check 4: "
        << (std::distance(FI(std::bind(predicate, std::placeholders::_1, 3))
                            .set_to_next_positive(tria.begin_active()),
                          FI(std::bind(predicate, std::placeholders::_1, 3),
                             tria.end())) ==
                static_cast<signed int>(tria.n_active_cells(3)) ?
              "OK" :
              "Failed")
        << std::endl;
    };


  // check 5: check that we loop over
  // all cells with a given subdomain
  // id
  if (true)
    {
      using FI = FilteredIterator<active_cell_iterator>;
      const IteratorFilters::SubdomainEqualTo predicate(1);
      FI                                      cell(predicate);
      cell.set_to_next_positive(tria.begin_active());
      active_cell_iterator endc(tria.end());
      active_cell_iterator cell1 = tria.begin_active();

      while (cell1->subdomain_id() != 1)
        ++cell1;

      while (true)
        {
          // move filtered iterator ahead
          ++cell;
          // move unfiltered iterator
          // ahead
          ++cell1;
          while ((cell1 != endc) && (cell1->subdomain_id() != 1))
            ++cell1;

          AssertThrow(cell == cell1, ExcInternalError());
          AssertThrow(cell1 == cell, ExcInternalError());

          if (cell.state() != IteratorState::valid)
            break;
        };
      AssertThrow(cell == endc, ExcInternalError());
      AssertThrow(cell1 == endc, ExcInternalError());

      logfile << "Check 5: OK" << std::endl;
    };
}


int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(4);
  test();

  return 0;
}
