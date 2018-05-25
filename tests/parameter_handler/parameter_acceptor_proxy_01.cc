//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2018 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal.II distribution.
//
//-----------------------------------------------------------



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
