// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Common subexpression elimination: Very simple example to demonstrate
// how the cascading levels of elimination would be created

#include <symengine/basic.h>
#include <symengine/parser.h>
#include <symengine/real_double.h>
#include <symengine/symbol.h>
#include <symengine/symengine_rcp.h>

#include <algorithm>
#include <map>
#include <sstream>
#include <vector>

#include "../tests.h"

using namespace dealii;
namespace SE = SymEngine;

int
main(int argc, char *argv[])
{
  initlog();

  SE::RCP<const SE::Symbol> a = SE::symbol("a");
  SE::RCP<const SE::Symbol> b = SE::symbol("b");
  SE::RCP<const SE::Symbol> c = SE::symbol("c");

  SE::vec_basic independents = {a, b, c};

  SE::vec_basic dependents = {SE::parse("a+b"),
                              SE::parse("(a+b)+(a+b)"),
                              SE::parse("(a+b)*c"),
                              SE::parse("(a+b)*(a+b)"),
                              SE::parse("pow(a+b,2)"),
                              SE::parse("pow(2,a+b)"),
                              SE::parse("pow(a+b,a+b)"),
                              SE::parse("pow(a+b,pow(a+b,a+b))"),
                              SE::parse("exp(a+b)"),
                              SE::parse("log((a+b)*c)"),
                              SE::parse("sin((a+b)*c)"),
                              SE::parse("asin((a+b)*c)")};

  deallog.push("Independents");
  for (unsigned int i = 0; i < independents.size(); i++)
    deallog << *(independents[i]) << std::endl;
  deallog.pop();
  deallog.push("Dependents");
  for (unsigned int i = 0; i < dependents.size(); i++)
    deallog << *(dependents[i]) << std::endl;
  deallog.pop();

  deallog.push("CSE elimination");
  SE::vec_pair  intermediate_symbols_exprs;
  SE::vec_basic reduced_exprs;
  SE::cse(intermediate_symbols_exprs, reduced_exprs, dependents);
  deallog.push("Intermediate reduced expressions");
  for (unsigned i = 0; i < intermediate_symbols_exprs.size(); ++i)
    {
      const SE::RCP<const SE::Basic> &cse_symbol =
        intermediate_symbols_exprs[i].first;
      const SE::RCP<const SE::Basic> &cse_expr =
        intermediate_symbols_exprs[i].second;
      deallog << i << ": " << *cse_symbol << " = " << *cse_expr << std::endl;
    }
  deallog.pop();
  deallog.push("Final reduced expressions for dependent variables");
  for (unsigned i = 0; i < reduced_exprs.size(); ++i)
    deallog << i << ": " << *(reduced_exprs[i]) << std::endl;
  deallog.pop();

  deallog << "OK" << std::endl;

  return 0;
}
