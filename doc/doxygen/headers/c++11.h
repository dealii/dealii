// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2017 by the deal.II authors
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


/**
 * @defgroup CPP11 deal.II and the C++11 standard
 *
 * Since version 9.0, deal.II requires a compiler that supports at
 * least <a href="http://en.wikipedia.org/wiki/C%2B%2B11">C++11</a>.
 * As part of this, many places in the internal implementation of
 * deal.II are now using features that were only introduced in C++11.
 * That said, deal.II also has functions and classes that make using
 * it with C++11 features easier.
 * 
 * One example is support for C++11
 * <a href="http://en.wikipedia.org/wiki/C++11#Range-based_for_loop">range-based
 * for loops</a>. deal.II-based codes often have many loops of the kind
 * @code
 *   Triangulation<dim> triangulation;
 *   ...
 *   typename Triangulation<dim>::active_cell_iterator
 *     cell = triangulation.begin_active(),
 *     endc = triangulation.end();
 *   for (; cell!=endc; ++cell)
 *     cell->set_refine_flag();
 * @endcode
 * Using C++11's range-based for loops, you can now write this as follows:
 * @code
 *   Triangulation<dim> triangulation;
 *   ...
 *   for (auto cell : triangulation.active_cell_iterators())
 *     cell->set_refine_flag();
 * @endcode
 * This relies on functions such as Triangulation::active_cell_iterators(),
 * and equivalents in the DoF handler classes,
 * DoFHandler::active_cell_iterators(), hp::DoFHandler::active_cell_iterators().
 * There are variants of these functions that provide iterator ranges
 * for all cells (not just the active ones) and for cells on individual
 * levels.
 */
