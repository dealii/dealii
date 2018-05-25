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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
