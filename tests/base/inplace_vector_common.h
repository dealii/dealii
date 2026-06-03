// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_tests_inplace_vector_common_h
#define dealii_tests_inplace_vector_common_h

#include <deal.II/base/std_cxx26/inplace_vector.h>

#include <ostream>

#include "../tests.h"

/**
 * Helper class for tracking object creation and destruction in inplace_vector.
 */
struct A
{
  static bool &
  logging()
  {
    static bool logging = true;
    return logging;
  }

  static int &
  n_ctors()
  {
    static int n_constructor_calls = 0;
    return n_constructor_calls;
  }

  static int &
  n_dtors()
  {
    static int n_destructor_calls = 0;
    return n_destructor_calls;
  }

  A()
  {
    ++n_ctors();
    if (logging())
      deallog << "A::A()" << std::endl;
  }

  A(A &&)
  {
    ++n_ctors();
    if (logging())
      deallog << "A::A(A&&)" << std::endl;
  }

  A(const A &)
  {
    ++n_ctors();
    if (logging())
      deallog << "A::A(const A&)" << std::endl;
  }

  ~A()
  {
    ++n_dtors();
    if (logging())
      deallog << "A::~A()" << std::endl;
  }

  A &
  operator=(const A &)
  {
    if (logging())
      deallog << "A::operator=(const A&)" << std::endl;
    return *this;
  }

  A &
  operator=(A &&)
  {
    if (logging())
      deallog << "A::operator=(A&&)" << std::endl;
    return *this;
  }

  bool
  operator==(const A &) const noexcept
  {
    return true;
  }

  bool
  operator!=(const A &) const noexcept
  {
    return false;
  }

  bool
  operator<(const A &) const noexcept
  {
    return false;
  }

  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int)
  {}
};

std::ostream &
operator<<(std::ostream &out, const A &)
{
  out << "{}";
  return out;
}

inline void
swap(A &, A &) noexcept
{
  // swapping is a no-op
}

template <typename T>
void
print_counts()
{
  if constexpr (std::is_same_v<T, A>)
    deallog << "ctors = " << A::n_ctors() << " dtors = " << A::n_dtors()
            << std::endl;
}

template <typename T>
void
print(const T &t)
{
  deallog << t;
}

template <typename T>
void
print(const std::vector<T> &vec)
{
  deallog << "{";
  if (vec.size() > 0)
    {
      print(vec[0]);
      for (std::size_t i = 1; i < vec.size(); ++i)
        {
          deallog << ", ";
          print(vec[i]);
        }
    }
  deallog << "}";
}

template <typename T, std::size_t N>
void
print(const std_cxx26::inplace_vector<T, N> &vec)
{
  if constexpr (!std::is_same_v<T, A>)
    {
      if (vec.size() > 0)
        {
          print(vec[0]);
          for (std::size_t i = 1; i < vec.size(); ++i)
            {
              deallog << ", ";
              print(vec[i]);
            }
        }
    }
  else
    deallog << "size : " << vec.size();
}

#endif
