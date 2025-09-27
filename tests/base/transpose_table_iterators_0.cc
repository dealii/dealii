// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check some basic behavior of the TransposeTable iterators
#include <deal.II/base/table.h>

#include <boost/core/demangle.hpp>

#include "../tests.h"

template <typename Container>
void
test()
{
  Container table(3, 4);

  auto it = table.begin();
  auto t2 = it;
  deallog << std::boolalpha;
  deallog << boost::core::demangle(typeid(Container).name()) << std::endl;
  deallog << "++it == it:           " << (++table.begin() == table.begin())
          << std::endl;
  deallog << "++it != it:           " << (++table.begin() != table.begin())
          << std::endl;
  deallog << "it++ == it:           " << (table.begin()++ == table.begin())
          << std::endl;
  deallog << "it++ != it + 1:       " << (table.begin()++ != table.begin() + 1)
          << std::endl;
  deallog << "++(++it) == it + 2:   "
          << (++(++table.begin()) == table.begin() + 2) << std::endl;
  deallog << "it++ == it + 1:       " << (table.begin()++ == table.begin() + 1)
          << std::endl;
  deallog << "--(++it) == it:       " << (--(++table.begin()) == table.begin())
          << std::endl;
  deallog << "(it--)-- == it:       " << ((table.end()--)-- == table.end())
          << std::endl;
  deallog << "end - 5 <= end:       " << ((table.end() - 5) <= table.end())
          << std::endl;
  deallog << "end <= end:           " << (table.end() <= table.end())
          << std::endl;
  deallog << "end >= end:           " << (table.end() >= table.end())
          << std::endl;
  deallog << "begin + 5 >= begin:   " << ((table.begin() + 5) >= table.begin())
          << std::endl;
  deallog << "begin + 5 <= begin:   " << ((table.begin() + 5) <= table.begin())
          << std::endl;
  deallog << "end - begin:          " << table.end() - table.begin()
          << std::endl;
  deallog << "begin - end:          " << table.begin() - table.end()
          << std::endl;
  deallog << "(begin + 2) - begin:  " << (table.begin() + 2) - table.begin()
          << std::endl;
  deallog << "begin + 6 < begin:    " << ((table.begin() + 6) < table.begin())
          << std::endl;
  deallog << "begin < begin + 1:    " << (table.begin() < (table.begin() + 6))
          << std::endl;
  deallog << "end - 5 < begin:      " << ((table.end() - 6) < table.begin())
          << std::endl;
  t2 = it + 5;
  deallog << "it+5 == (t2 = (it+5)):" << (it + 5 == t2) << std::endl;
  const auto it2 = table.end() - 5;
  deallog << "end - 5 position:     " << it2->row() << ", " << it2->column()
          << std::endl;
  const auto it3 = table.begin() + 5;
  deallog << "begin + 5 position:   " << it3->row() << ", " << it3->column()
          << std::endl;
}

int
main()
{
  initlog();

  // Test an empty table
  {
    TransposeTable<double> transpose_table;
    Assert(transpose_table.begin() == transpose_table.end(),
           ExcMessage("The beginning and end iterators should be equal for an "
                      "empty table."));
    TransposeTable<double>::const_iterator begin = transpose_table.begin();
    TransposeTable<double>::const_iterator end   = transpose_table.end();
    Assert(begin == end,
           ExcMessage("The beginning and end const iterators should"
                      " be equal for an empty table."));
  }

  // test some things with accessors
  {
    TransposeTable<double> table(2, 2);

    MatrixTableIterators::Accessor<TransposeTable<double>,
                                   true,
                                   MatrixTableIterators::Storage::column_major>
      a3(&table, 2);
    MatrixTableIterators::Accessor<TransposeTable<double>,
                                   false,
                                   MatrixTableIterators::Storage::column_major>
      a4(&table, 2);
    deallog << "Accessors refer to the same entry: "
            << (&(MatrixTableIterators::Accessor<
                    TransposeTable<double>,
                    true,
                    MatrixTableIterators::Storage::column_major>(a4)
                    .value()) == &(a3.value()))
            << std::endl;
  }

  // test a non-empty rectangular table
  test<TransposeTable<double>>();
  test<const TransposeTable<TransposeTable<double>>>();

  {
    // check that we can convert
    TransposeTable<double>                 table(3, 4);
    TransposeTable<double>::iterator       it0 = table.begin();
    TransposeTable<double>::const_iterator it1 = it0;
    deallog << "converted iterators are equal: " << (it1 == it0) << std::endl;
  }

  deallog << "OK" << std::endl;
}
