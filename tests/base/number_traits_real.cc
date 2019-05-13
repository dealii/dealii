// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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


// check numbers::NumberTraits for real data types

#include <limits>
#include <typeinfo>

#include "../tests.h"


template <typename number>
void
check(const number &x)
{
  deallog << "typeid(x).name() = " << typeid(x).name() << std::endl;

  deallog << "typeid(NumberTraits<number>::real_type).name() = "
          << typeid(typename numbers::NumberTraits<number>::real_type).name()
          << std::endl;

  deallog << numbers::NumberTraits<number>::conjugate(x) << std::endl;

  deallog << numbers::NumberTraits<number>::abs_square(x) << std::endl;

  deallog << numbers::NumberTraits<number>::abs(x) << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check((float)1.5);
  check((float)-1.5);

  check((double)1.5);
  check((double)-1.5);

  check((long double)1.5);
  check((long double)-1.5);

  return 0;
}
