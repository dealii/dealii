// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * @defgroup CPP11 deal.II and Modern C++ standards
 *
 * Since version 9.6, deal.II requires a compiler that supports at least
 * <a href="https://en.wikipedia.org/wiki/C%2B%2B17">C++17</a>. Large parts
 * of the library now depend on modern language constructs which are
 * documented here.
 *
 * One example is support for C++11
 * <a href="https://en.wikipedia.org/wiki/C++11#Range-based_for_loop">range-based
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
 * This works in the same way with Triangulation::active_cell_iterators()
 * and DoFHandler::active_cell_iterators().
 * There are variants of these functions that provide iterator ranges
 * for all cells (not just the active ones) and for cells on individual
 * levels.
 *
 * There are numerous other functions in the library that allow for
 * the idiomatic use of range-based for loops. Examples are
 * ReferenceCell::face_indices(), ReferenceCell::vertex_indices(),
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
 * Some functions such as determinant() are specified as `constexpr`: these rely
 * on the generalized constexpr support available in C++14. Some functions,
 * such as unit_symmetric_tensor(), rely on further developments of `constexpr`
 * only available in C++17 and newer. As such, this function is declared as
 * @code
 * template <int dim, typename Number>
 * DEAL_II_CONSTEXPR inline SymmetricTensor<2, dim, Number>
 * unit_symmetric_tensor();
 * @endcode
 * The macro @ref DEAL_II_CONSTEXPR expands to `constexpr` if the compiler
 * supports enough `constexpr` features (such as loops). If the compiler does
 * not then this macro expands to nothing.
 *
 * Functions declared as `constexpr` can be evaluated at compile time. Hence code
 * like
 * @code
 * constexpr double det_A = determinant(A);
 * @endcode
 * assuming `A` is declared with the `constexpr` specifier, will typically
 * result in compile-time constants. This example shows the performance gains of
 * using `constexpr` because here we performed an operation with
 * $O(\text{dim}^3)$ complexity during compile time, avoiding any runtime cost.
 */



/**
 * Similarly, deal.II defined the C++17 library features it used before
 * requiring C++17 in this namespac.
 *
 * The most notable example is the <a
 * href="https://en.cppreference.com/w/cpp/utility/optional">`std::optional`</a> class
 * that was introduced to C++ starting with the C++17 standard.
 */
namespace std_cxx17
{}



/**
 * deal.II currently only requires a C++17-conforming compiler, but there are a
 * number of functions and classes from the C++20 standard that are easy to
 * provide also in case the compiler only supports C++17. These are collected
 * in the current namespace.
 *
 * One example is the <a
 * href="https://en.cppreference.com/w/cpp/ranges/iota_view">`std::ranges::iota_view`</a>
 * class that was introduced to C++ starting with the C++20
 * standard. It is used as the return type for the
 * ReferenceCell::face_indices(), ReferenceCell::vertex_indices(), and
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
