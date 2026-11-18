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


// Test AlignedVector::insert() when called repeatedly: it should reuse the
// existing capacity and only reallocate a logarithmic number of times, rather
// than reallocating (and moving the whole array) on every call. Also check
// correctness against std::vector for trivially and non-trivially copyable
// types, for different insertion positions, and for different iterator types.

#include <deal.II/base/aligned_vector.h>

#include <string>
#include <vector>

#include "../tests.h"


// A non-trivially copyable type that tracks the number of live objects so
// that we can detect leaks and double destructions.
struct Tracked
{
  static int n_alive;

  Tracked(const int v = 0)
    : value(std::to_string(v))
  {
    ++n_alive;
  }

  Tracked(const Tracked &other)
    : value(other.value)
  {
    ++n_alive;
  }

  Tracked(Tracked &&other) noexcept
    : value(std::move(other.value))
  {
    ++n_alive;
  }

  Tracked &
  operator=(const Tracked &) = default;

  ~Tracked()
  {
    --n_alive;
  }

  bool
  operator==(const Tracked &other) const
  {
    return value == other.value;
  }

  std::string value;
};

int Tracked::n_alive = 0;



template <typename T>
void
check_equal(const AlignedVector<T> &a, const std::vector<T> &b)
{
  AssertThrow(a.size() == b.size(), ExcInternalError());
  for (std::size_t i = 0; i < b.size(); ++i)
    AssertThrow(a[i] == b[i], ExcInternalError());
}



template <typename T>
void
test_repeated_insert()
{
  const unsigned int n = 10000;

  // insert at the end, one element at a time
  {
    AlignedVector<T> a;
    std::vector<T>   b;
    unsigned int     n_reallocations = 0;
    for (unsigned int i = 0; i < n; ++i)
      {
        const T    value(i);
        const auto old_data = a.data();
        const auto it       = a.insert(a.end(), &value, &value + 1);
        b.insert(b.end(), value);
        AssertThrow(it == a.end() - 1, ExcInternalError());
        if (a.data() != old_data)
          ++n_reallocations;
      }
    check_equal(a, b);
    deallog << "insert at end: size " << a.size()
            << ", few reallocations: " << (n_reallocations <= 20) << std::endl;
  }

  // insert at the front and in the middle, one and several elements at a
  // time
  {
    AlignedVector<T> a;
    std::vector<T>   b;
    for (unsigned int i = 0; i < 500; ++i)
      {
        const T    values[3] = {T(i), T(i + 1000), T(i + 2000)};
        const auto n_values  = 1 + i % 3;

        auto it = a.insert(a.begin(), values, values + n_values);
        b.insert(b.begin(), values, values + n_values);
        AssertThrow(it == a.begin(), ExcInternalError());

        const auto offset = a.size() / 2;
        it = a.insert(a.begin() + offset, values, values + n_values);
        b.insert(b.begin() + offset, values, values + n_values);
        AssertThrow(it == a.begin() + offset, ExcInternalError());
      }
    check_equal(a, b);
    deallog << "insert at front and middle: size " << a.size() << std::endl;
  }

  // insert ranges from other containers: pointers into another AlignedVector
  // and iterators of a std::vector
  {
    AlignedVector<T> a;
    std::vector<T>   b;
    AlignedVector<T> other;
    std::vector<T>   other_std;
    for (unsigned int i = 0; i < 8; ++i)
      {
        other.push_back(T(100 + i));
        other_std.push_back(T(200 + i));
      }
    for (unsigned int i = 0; i < 50; ++i)
      {
        const auto offset = a.size() / 3;
        a.insert(a.begin() + offset, other.begin(), other.end());
        b.insert(b.begin() + offset, other.begin(), other.end());
        a.insert(a.begin() + offset, other_std.begin() + 2, other_std.end());
        b.insert(b.begin() + offset, other_std.begin() + 2, other_std.end());
      }
    check_equal(a, b);
    deallog << "insert from other containers: size " << a.size() << std::endl;
  }

  // inserting an empty range should do nothing
  {
    AlignedVector<T>     a(3);
    const std::vector<T> empty;
    const auto           old_data = a.data();
    const auto it = a.insert(a.begin() + 1, empty.begin(), empty.end());
    AssertThrow(it == a.begin() + 1, ExcInternalError());
    AssertThrow(a.size() == 3, ExcInternalError());
    AssertThrow(a.data() == old_data, ExcInternalError());
  }
}



int
main()
{
  initlog();

  deallog << "unsigned int" << std::endl;
  test_repeated_insert<unsigned int>();

  deallog << "Tracked" << std::endl;
  test_repeated_insert<Tracked>();
  deallog << "live objects: " << Tracked::n_alive << std::endl;
}
