// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2023 by the deal.II Authors
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

#ifndef dealii_test_setup_h
#define dealii_test_setup_h

#include <ginkgo/ginkgo.hpp>


template <typename T>
gko::remove_complex<T>
get_relative_error(const T a, const T b)
{
  return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
}

template <typename T>
constexpr gko::remove_complex<T>
  tol = 10 * std::numeric_limits<gko::remove_complex<T>>::epsilon();


#define TEST_ASSERT(_assertion)                                     \
  if (!(_assertion))                                                \
    {                                                               \
      deallog << std::string(TEST_name) + ":FAIL with " #_assertion \
              << std::endl;                                         \
      self->test_success = false;                                   \
      return;                                                       \
    }                                                               \
  static_assert(true, "Enforce ; after macro")


#define TEST_ASSERT_NEAR(_a, _b, _tol)                                  \
  {                                                                     \
    auto rel_err = get_relative_error(_a, _b);                          \
    if (!(rel_err <= _tol))                                             \
      {                                                                 \
        deallog << std::string(TEST_name) + ":FAIL with rel_error(" #_a \
                                            ", " #_b ")["               \
                << rel_err << "] <= " << _tol << std::endl;             \
        self->test_success = false;                                     \
        return;                                                         \
      }                                                                 \
  }                                                                     \
  static_assert(true, "Enforce ; after macro")


#define TEST_ASSERT_THROW(_cmd, _expected_exc) \
  try                                          \
    {                                          \
      _cmd;                                    \
      self->test_success = false;              \
      return;                                  \
    }                                          \
  catch (const _expected_exc &)                \
    {}                                         \
  static_assert(true, "Enforce ; after macro")



#define TEST(_name)                                             \
  struct _name                                                  \
  {                                                             \
    _name()                                                     \
      : test_success(true)                                      \
    {                                                           \
      _name::run(this);                                         \
      if (test_success)                                         \
        deallog << std::string(TEST_name) + ":OK" << std::endl; \
    }                                                           \
                                                                \
    static void                                                 \
    run(_name *self);                                           \
                                                                \
    static constexpr char TEST_name[] = #_name;                 \
                                                                \
    bool test_success;                                          \
  };                                                            \
  constexpr char _name::TEST_name[];                            \
  void           _name::run(_name *self)


#endif // dealii_test_setup_h
