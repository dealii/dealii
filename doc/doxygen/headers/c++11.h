// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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


/**
 * @defgroup CPP11 deal.II and the C++11 standard
 *
 * At present, deal.II only requires a compiler that conforms to the
 * <a href="http://en.wikipedia.org/wiki/C%2B%2B#Standardization">C++98</a>
 * standard and does not rely on compilers to either
 * provide the features introduced in
 * <a href="http://en.wikipedia.org/wiki/C%2B%2B03">C++03</a> or
 * <a href="http://en.wikipedia.org/wiki/C%2B%2B11">C++11</a>
 * 
 * That said, deal.II interfaces with C++11 in several ways as
 * outlined below.
 * 
 * 
 * <h3>Use of C++11 classes and substitution by BOOST</h3>
 *
 * deal.II makes use of many of the classes that were only
 * added as part of C++11. This includes std::shared_ptr,
 * std::function, std::bind, std::tuple and a number of others.
 * Because we do not assume that the compiler actually supports
 * C++11, there needs to be a way to ensure that these classes
 * are available also for pre-C++11 compilers. This is done using
 * the following approach:
 *
 * - We create a namespace std_cxx11.
 * - If the compiler supports C++11, we import the relevant classes
 *   and functions into this namespace using statements such as
 *   @code
 *     namespace std_cxx11 {  using std::shared_ptr;  }
 *   @endcode
 * - If the compiler does not support C++11, if its support for
 *   C++11 is incomplete, or if it is buggy, then we use as a fallback
 *   the corresponding classes and functions provided by the
 *   <a href="http://www.boost.org">BOOST library</a> through
 *   statements such as
 *   @code
 *     namespace std_cxx11 {  using boost::shared_ptr;  }
 *   @endcode
 *
 * Consequently, namespace std_cxx11 contains all of the symbols
 * we require. The classes that can be used this way are obviously
 * a subset of the intersection between C++11 and what BOOST provides.
 *
 *
 * <h3>Support for C++11 range-based for loops</h3>
 *
 * C++11 provides many new core language features, such as
 * rvalue references and move semantics, initialized lists, tuples,
 * variadic templates and
 * others. For a complete list, see  http://en.wikipedia.org/wiki/C++11 .
 * We can not use most of these in deal.II itself because we cannot rely
 * on compilers supporting them.
 *
 * However, this does not preclude users from using such features in their
 * own applications if they can be reasonably sure that the compilers on
 * all of the systems they will work on do support C++11. An example are
 * <a href="http://en.wikipedia.org/wiki/C++11#Type_inference">automatically
 * typed variables</a>.
 *
 * deal.II does provide some features that make programming simpler when using
 * C++11. This is true, in particular, for
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
