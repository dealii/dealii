// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#ifndef dealii_cxx20_type_traits_h
#define dealii_cxx20_type_traits_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_HAVE_CXX20
#  include <type_traits>
#endif

DEAL_II_NAMESPACE_OPEN

namespace std_cxx20
{
#ifdef DEAL_II_HAVE_CXX20
  using std::type_identity;
  using std::type_identity_t;
#else
  /**
   * A template class that simply exports its template argument as a local
   * alias. This class, while at first appearing useless, makes sense in the
   * following context: if you have a function template as follows:
   * @code
   * template <typename T>
   * void f(T, T);
   * @endcode
   * then it can't be called in an expression like <code>f(1, 3.141)</code>
   * because the type <code>T</code> of the template can not be deduced in a
   * unique way from the types of the arguments. However, if the template is
   * written as
   * @code
   * template <typename T>
   * void f(T, std_cxx20::type_identity_t<T>);
   * @endcode
   * then the call becomes valid: the type <code>T</code> is not deducible
   * from the second argument to the function, so only the first argument
   * participates in template type resolution.
   *
   * The context for this feature is as follows: consider
   * @code
   * template <typename RT, typename A>
   * void forward_call(RT (*p) (A), A a)
   * {
   *   p(a);
   * }
   *
   * void h (double);
   *
   * void g()
   * {
   *   forward_call(&h, 1);
   * }
   * @endcode
   * This code fails to compile because the compiler can't decide whether the
   * template type <code>A</code> should be <code>double</code> (from the
   * signature of the function given as first argument to
   * <code>forward_call</code>, or <code>int</code> because the expression
   * <code>1</code> has that type. Of course, what we would like the compiler
   * to do is simply cast the <code>1</code> to <code>double</code>. We can
   * achieve this by writing the code as follows:
   * @code
   * template <typename RT, typename A>
   * void forward_call(RT (*p) (A), std_cxx20::type_identity_t<A> a)
   * {
   *   p(a);
   * }
   *
   * void h (double);
   *
   * void g()
   * {
   *   forward_call(&h, 1);
   * }
   * @endcode
   */
  template <typename T>
  struct type_identity
  {
    using type = T;
  };

  template <class T>
  using type_identity_t = typename type_identity<T>::type;

#endif
} // namespace std_cxx20

DEAL_II_NAMESPACE_CLOSE

#endif
