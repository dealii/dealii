// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test DynamicSparsityPattern::iterator and document bug with std::vector
// iterators

/*
crash:
/usr/include/c++/4.8/debug/safe_iterator.h:506:error: attempt to compare a
    singular iterator to a singular iterator.

Objects involved in the operation:
iterator "lhs" @ 0x0x7fffffffc860 {
type =
N11__gnu_debug14_Safe_iteratorIN9__gnu_cxx17__normal_iteratorIPKjSt6vectorIjSaIjEEEENSt7__debug6vectorIjS6_EEEE
(constant iterator); state = singular;
}
iterator "rhs" @ 0x0x7fffffffd670 {
type =
N11__gnu_debug14_Safe_iteratorIN9__gnu_cxx17__normal_iteratorIPKjSt6vectorIjSaIjEEEENSt7__debug6vectorIjS6_EEEE
(constant iterator); state = singular;
}

backtrace:
#3  0x0000000000416c29 in
__gnu_debug::operator==<__gnu_cxx::__normal_iterator<unsigned int const*,
std::vector<unsigned int, std::allocator<unsigned int> > >,
std::__debug::vector<unsigned int, std::allocator<unsigned int> > >
(__lhs=<error reading variable: Cannot access memory at address 0x0>,
__rhs=<error reading variable: Cannot access memory at address 0x0>) at
/usr/include/c++/4.8/debug/safe_iterator.h:503 #4  0x000000000041075b in
operator== (other=..., this=0x7fffffffc850) at
/scratch/deal-git/include/deal.II/lac/dynamic_sparsity_pattern.h:743 #5
operator== (other=..., this=0x7fffffffc850) at
/scratch/deal-git/include/deal.II/lac/dynamic_sparsity_pattern.h:868 #6
operator!= (other=..., this=0x7fffffffc850) at
/scratch/deal-git/include/deal.II/lac/dynamic_sparsity_pattern.h:877 #7  iterate
(sp=...) at
/scratch/deal-git/tests/bits/dynamic_sparsity_pattern_iterator_04.cc:30
*/

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
iterate(DynamicSparsityPattern &sp)
{
  DynamicSparsityPattern::const_iterator i = sp.begin();
  for (; i != sp.end(); ++i)
    deallog << i->row() << ' ' << i->column() << std::endl;

  deallog << "OK" << std::endl;

  {
    for (unsigned int row = 0; row < sp.n_rows(); ++row)
      {
        DynamicSparsityPattern::iterator col     = sp.begin(row),
                                         end_col = sp.end(row);
        deallog << "row " << row << ':' << std::endl;
        for (; col != end_col; ++col)
          {
            deallog << "row= " << col->row() << ", col= " << col->column()
                    << ", index= " << col->index() << std::endl;
          }
      }
  }
  deallog << "OK" << std::endl;
}


void
test()
{
  {
    DynamicSparsityPattern sp(5, 4);
    sp.add(0, 0);
    sp.add(0, 1);
    sp.add(3, 3);
    sp.compress();

    deallog << "** 5 by 4 ** " << std::endl;
    iterate(sp);
  }
  {
    DynamicSparsityPattern sp(5, 5);
    sp.add(0, 0);
    sp.add(0, 1);
    sp.add(3, 3);
    sp.compress();

    deallog << "** 5 by 5 ** " << std::endl;
    iterate(sp);
  }
}



int
main()
{
  initlog();

  test();
}
