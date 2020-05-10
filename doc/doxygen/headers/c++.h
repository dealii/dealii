// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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
 *   for (auto &cell : triangulation.active_cell_iterators())
 *     cell->set_refine_flag();
 * @endcode
 * This relies on functions such as Triangulation::active_cell_iterators(),
 * and equivalents in the DoF handler classes,
 * DoFHandler::active_cell_iterators(), hp::DoFHandler::active_cell_iterators().
 * There are variants of these functions that provide iterator ranges
 * for all cells (not just the active ones) and for cells on individual
 * levels.
 *
 * There are numerous other functions in the library that allow for
 * the idiomatic use of range-based for loops. Examples are
 * GeometryInfo::face_indices(), GeometryInfo::vertex_indices(),
 * FEValuesBase::quadrature_point_indices(), among many others.
 *
 * C++11 also introduces the concept of
 * [constexpr](https://en.cppreference.com/w/cpp/language/constexpr)
 * variables and function. The variables defined as `constexpr` are constant
 * values that are computed during the compilation of the program and therefore
 * have zero runtime cost associated with their initialization. Additionally,
 * `constexpr` constants have properly defined lifetimes which prevent the
 * so-called "static initialization order fiasco" completely. %Functions can be
 * marked as `constexpr`, indicating that they can produce compile-time
 * constant return values if their input arguments are constant expressions.
 * Additionally, classes with at least one `constexpr` constructor can be
 * initialized as `constexpr`.
 *
 * As an example, since the constructor Tensor::Tensor(const array_type &) is
 * `constexpr`, we can initialize a tensor with an array during compile time
 * as:
 * @code
 * constexpr double[2][2] entries = {{1., 0.}, {0., 1.}};
 * constexpr Tensor<2, 2> A(entries);
 * @endcode
 * Here, the contents of A are not stored on the stack. Rather, they are
 * initialized during compile time and inserted into the `.data` portion
 * of the executable program. The program can use these values at runtime
 * without spending time for initialization. Initializing tensors can be
 * simplified in one line.
 * @code
 * constexpr Tensor<2, 2> A({{1., 0.}, {0., 1.}});
 * @endcode
 * Some functions such as determinant() are specified as `constexpr` but they
 * require a compiler with C++14 capability. As such, this function is
 * internally declared as:
 * @code
 * template <int dim, typename Number>
 * DEAL_II_CONSTEXPR Number determinant(const Tensor<2, dim, Number> &t);
 * @endcode
 * The macro @ref DEAL_II_CONSTEXPR simplifies to `constexpr` if a C++14-capable
 * compiler is available. Otherwise, for old compilers, it ignores
 * DEAL_II_CONSTEXPR altogether.
 * Therefore, with newer compilers, the user can write
 * @code
 * constexpr double det_A = determinant(A);
 * @endcode
 * assuming `A` is declared with the `constexpr` specifier. This example shows
 * the performance gains of using `constexpr` because here we performed an
 * operation with $O(\text{dim}^3)$ complexity during compile time, avoiding
 * any runtime cost.
 */



/**
 * deal.II currently only requires a C++11-conforming compiler, but there are a
 * number of functions and classes from the C++14 standard that are easy to
 * provide also in case the compiler only supports C++11. These are collected
 * in the current namespace.
 *
 * The most notable example is the <a
 * href="https://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique">`std::make_unique`</a>
 * function which is arguably an oversight for not having been
 * included in C++11 (given that there is <a
 * href="https://en.cppreference.com/w/cpp/memory/shared_ptr/make_shared">`std::make_shared`</a>
 * in C++11).
 *
 * There are other small additions in this namespace that allow us to
 * use C++14 features at this point already, even though we don't
 * require a C++14-compliant compiler.
 *
 * @note If the compiler in use actually does support C++14, then the
 *   contents of this namespace are simply imported classes and
 *   functions from namespace `std`. That is, we fall back to what the
 *   compiler provides, rather than our own implementations.
 */
namespace std_cxx14
{}



/**
 * deal.II currently only requires a C++11-conforming compiler, but there are a
 * number of functions and classes from the C++17 standard that are easy to
 * provide also in case the compiler only supports C++11. These are collected
 * in the current namespace.
 *
 * The most notable example is the <a
 * href="https://en.cppreference.com/w/cpp/utility/optional">`std::optional`</a> class
 * that was introduced to C++ starting with the C++17 standard.
 *
 * There are other small additions in this namespace that allow us to
 * use C++17 features at this point already, even though we don't
 * require a C++17-compliant compiler.
 *
 * @note If the compiler in use actually does support C++17, then the
 *   contents of this namespace are simply imported classes and
 *   functions from namespace `std`. That is, we fall back to what the
 *   compiler provides, rather than our own implementations.
 */
namespace std_cxx17
{}



/**
 * deal.II currently only requires a C++11-conforming compiler, but there are a
 * number of functions and classes from the C++20 standard that are easy to
 * provide also in case the compiler only supports C++11. These are collected
 * in the current namespace.
 *
 * One example is the <a
 * href="https://en.cppreference.com/w/cpp/ranges/iota_view">`std::ranges::iota_view`</a>
 * class that was introduced to C++ starting with the C++20
 * standard. It is used as the return type for the
 * GeometryInfo::face_indices(), GeometryInfo::vertex_indices(), and
 * FEValuesBase::quadrature_point_indices() functions, among others,
 * to support range-based for loops (see @ref CPP11 for examples of
 * range-based for loops, as well as the documentation of the
 * functions mentioned above).
 *
 * There are other small additions in this namespace that allow us to
 * use C++20 features at this point already, even though we don't
 * require a C++20-compliant compiler.
 *
 * @note If the compiler in use actually does support C++20, then the
 *   contents of this namespace are simply imported classes and
 *   functions from namespace `std`. That is, we fall back to what the
 *   compiler provides, rather than our own implementations.
 */
namespace std_cxx20
{}
