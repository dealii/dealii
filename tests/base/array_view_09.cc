// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


// test the make_array_view functions that create ArrayViews from iterators

#include <deal.II/base/array_view.h>

#include <boost/container/static_vector.hpp>

#include <type_traits>

#include "../tests.h"

template <typename T>
constexpr bool
is_const_reference()
{
  return std::is_reference<T>::value &&
         std::is_const<typename std::remove_reference<T>::type>::value;
}

void
test()
{
  // test on a built-in array
  {
    int  v[10];
    auto a = make_array_view(std::begin(v) + 1, std::end(v));
    AssertThrow(a.begin() == v + 1, ExcInternalError());
    AssertThrow(a.end() == v + sizeof(v) / sizeof(v[0]), ExcInternalError());
    AssertThrow(a.begin() + 2 == &v[3], ExcInternalError());
    a[2] = 42;
    AssertThrow(v[3] == 42, ExcInternalError());

    // check that we cannot create a backwards array
    try
      {
        make_array_view(std::end(v), std::begin(v));
      }
    catch (const ExceptionBase &exc)
      {
        deallog << exc.get_exc_name() << std::endl;
      }
  }

  // test on the std::array class template
  {
    std::array<double, 10> v;
    std::fill(v.begin(), v.end(), 42.0);
    const std::array<double, 10> &v2 = v;
    const auto                    a  = make_array_view(v2.begin(), v2.end());
    static_assert(is_const_reference<decltype(*a.begin())>(),
                  "type should be const");
    static_assert(is_const_reference<decltype(*a.end())>(),
                  "type should be const");

    v[5] = 10;
    AssertThrow(a[5] == 10, ExcInternalError());
  }

  // test for a vector (preserve iterator constness)
  {
    std::vector<double> v(10);
    std::fill(v.begin(), v.end(), 42.0);
    auto a = make_array_view(v.cbegin() + 2, v.cend());
    // a should be ArrayView<const double>
    static_assert(!std::is_const<decltype(a)>::value,
                  "a should not be const (but has const value)");
    static_assert(std::is_const<decltype(a)::value_type>::value,
                  "a::value_type should be const");
    static_assert(is_const_reference<decltype(*a.begin())>(),
                  "type should be const");
    static_assert(is_const_reference<decltype(*a.end())>(),
                  "type should be const");
    v[2] = 10.0;
    AssertThrow(a[0] == v[2], ExcInternalError());
  }

  // test for a (new in boost 1.54) vector class
  {
    boost::container::static_vector<double, 20> v(10);
    std::fill(v.begin(), v.end(), 42.0);

    const auto a = make_array_view(v.cbegin() + 2, v.cend());
    AssertThrow(a.size() + 2 == v.size(), ExcInternalError());
    // some older versions of boost (1.57 in particular) do not implement
    // static_vector::cbegin() and static_vector::cend() correctly, so ignore
    // the type checking in that case
#if BOOST_VERSION >= 106200
    // a should be const ArrayView<const double>
    static_assert(std::is_const<decltype(a)>::value,
                  "a should not be const (but has const value)");
    static_assert(std::is_const<decltype(a)::value_type>::value,
                  "a::value_type should be const");
    static_assert(is_const_reference<decltype(*a.begin())>(),
                  "type should be const");
    static_assert(is_const_reference<decltype(*a.end())>(),
                  "type should be const");
#endif
    v[2] = 10.0;
    AssertThrow(a[0] == v[2], ExcInternalError());

    // check that we cannot create a backwards array
    try
      {
        make_array_view(std::end(v), std::begin(v));
      }
    catch (const ExceptionBase &exc)
      {
        deallog << exc.get_exc_name() << std::endl;
      }
      // we only trigger the exception in debug mode
      // just output the expected error message in release mode
      // to make CTest happy
#ifndef DEBUG
    deallog
      << "ExcMessage(\"The beginning of the array view should be before the end.\")"
      << std::endl;
    deallog
      << "ExcMessage(\"The beginning of the array view should be before the end.\")"
      << std::endl;
#endif
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  deal_II_exceptions::disable_abort_on_exception();
  test();
}
