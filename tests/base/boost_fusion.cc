// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// when updating BOOST to 1.56 in the summer of 2014, we thought we
// could get rid of BOOST.fusion. but this doesn't work for some users
// since boost/math/special_functions/erf.hpp uses it. this test is a
// reminder that we can't do this again in the future

#include <boost/math/special_functions/erf.hpp>

#include "../tests.h"


int
main()
{
  initlog();

  deallog << boost::math::erf(0.5) << std::endl;
}
