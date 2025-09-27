// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>

#include "../tests.h"

// Test proxy class


int
main()
{
  initlog();
  auto &prm = ParameterAcceptor::prm;

  ParameterAcceptorProxy<Functions::ParsedFunction<1>> f1("Function 1D");
  ParameterAcceptorProxy<Functions::ParsedFunction<2>> f2("Function 2D");
  ParameterAcceptorProxy<Functions::ParsedFunction<3>> f3("Function 3D");

  ParameterAcceptor::initialize(
    SOURCE_DIR
    "/parameter_acceptor_parameters/parameter_acceptor_proxy_01.prm");
  prm.log_parameters(deallog);
}
