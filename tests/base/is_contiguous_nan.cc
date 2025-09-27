// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that is_contiguous works correctly for arrays that contain
// signaling NaNs. In the first version of the is_contiguous() helper
// function, we compared the *objects* pointed to, as opposed to the
// *addresses* of these objects. This required (i) that these objects
// have an operator==, and (ii) that a==a for any object 'a'. The
// latter is not true if 'a' is a double that contains a
// signaling_nan<double>(), and consequently, the is_contiguous()
// function returned false for arrays with such entries.
//
// The failure of the function to recognize that the array elements
// are contiguous led to downstream failures of make_array_view().

#include <deal.II/base/array_view.h>
#include <deal.II/base/signaling_nan.h>

#include "../tests.h"

void
test()
{
  std::vector<double> tmp(2);
  tmp[0] = numbers::signaling_nan<double>();
  tmp[1] = numbers::signaling_nan<double>();
  deallog << (internal::ArrayViewHelper::is_contiguous(tmp.begin(), tmp.end()) ?
                "true" :
                "false")
          << std::endl;
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();
  initlog();

  test();
}
