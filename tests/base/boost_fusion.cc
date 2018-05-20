// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2017 by the deal.II authors
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

// when updating BOOST to 1.56 in the summer of 2014, we thought we
// could get rid of BOOST.fusion. but this doesn't work for some users
// since boost/math/special_functions/erf.hpp uses it. this test is a
// reminder that we can't do this again in the future

#include "../tests.h"

#include <boost/math/special_functions/erf.hpp>

int
main()
{
  initlog();

  deallog << boost::math::erf(0.5) << std::endl;
}
