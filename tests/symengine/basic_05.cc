// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
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

// Check that SymEngine can do some optimsation on diffentiated expressions
// using lambda functions
// This is an extension of basic_04.cc

// References:
// https://github.com/symengine/symengine/blob/master/symengine/tests/eval/test_lambda_double.cpp

#include "../tests.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra-semi"

#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/integer.h>
#include <symengine/lambda_double.h>
#include <symengine/llvm_double.h>
#include <symengine/mul.h>
#include <symengine/parser.h>
#include <symengine/pow.h>
#include <symengine/real_double.h>
#include <symengine/symbol.h>

#pragma GCC diagnostic pop

#include <chrono>
#include <iostream>
#include <vector>

namespace SE = SymEngine;

int
main(int argc, char *argv[])
{
  initlog();

  const unsigned int n_runs = 1000;

  SE::RCP<const SE::Symbol> t_11 = SE::symbol("t_11");
  SE::RCP<const SE::Symbol> t_12 = SE::symbol("t_12");
  SE::RCP<const SE::Symbol> s    = SE::symbol("s");
  SE::RCP<const SE::Symbol> t_02 = SE::symbol("t_02");
  SE::RCP<const SE::Symbol> v_0  = SE::symbol("v_0");
  SE::RCP<const SE::Symbol> v_1  = SE::symbol("v_1");
  SE::RCP<const SE::Symbol> v_2  = SE::symbol("v_2");
  SE::RCP<const SE::Symbol> t_10 = SE::symbol("t_10");
  SE::RCP<const SE::Symbol> t_00 = SE::symbol("t_00");
  SE::RCP<const SE::Symbol> t_01 = SE::symbol("t_01");
  SE::RCP<const SE::Symbol> t_20 = SE::symbol("t_20");
  SE::RCP<const SE::Symbol> t_21 = SE::symbol("t_21");
  SE::RCP<const SE::Symbol> t_22 = SE::symbol("t_22");

  SE::vec_basic v = {
    t_11, t_12, s, t_02, v_0, v_1, v_2, t_10, t_00, t_01, t_20, t_21, t_22};

  SE::RCP<const SE::Basic> h =
    SE::parse("2.*s**2.2*(v_0**2 + v_1**2 + v_2**2)**3*"
              "(-1.*t_20*(1.*t_01*t_12 - 1.*t_02*t_11) "
              "+ 1.*t_21*(1.*t_00*t_12 - 1.*t_02*t_10) "
              "- 1.*t_22*(1.*t_00*t_11 - 1.*t_01*t_10))"
              "*(1.*t_21*t_12 - 1.*t_22*t_11)");

  SE::vec_basic diffs = {h->diff(t_20),
                         h->diff(t_21),
                         h->diff(t_22),
                         h->diff(t_10),
                         h->diff(t_11),
                         h->diff(t_12),
                         h->diff(t_00),
                         h->diff(t_01),
                         h->diff(t_02),
                         h->diff(s)};

  std::cout << *h << std::endl;

  SE::map_basic_basic dict;
  std::vector<double> vals(v.size());
  for (unsigned i = 0; i < v.size(); i++)
    {
      dict[v[i]] = SE::real_double(i);
      vals[i]    = i;
    }

  double res[10] = {}, res1[10] = {}, res2[10] = {};
  double r = 0.0, r2 = 0.0, r3 = 0.0;

  // Standard substitution

  auto t1 = std::chrono::high_resolution_clock::now();
  for (unsigned j = 0; j < n_runs; j++)
    {
      for (unsigned k = 0; k < diffs.size(); k++)
        {
          res[k] = static_cast<const SE::RealDouble &>(*diffs[k]->subs(dict)).i;
          r += res[k];
        }
    }
  auto         t2 = std::chrono::high_resolution_clock::now();
  const double diff_symm_subs =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout
    << "Subs " << n_runs << " calls                  :"
    << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
    << " us" << std::endl;

  deallog << "Value (subs): " << r << std::endl;

  // Optimisation via lambda functions

  t1 = std::chrono::high_resolution_clock::now();
  SE::LambdaRealDoubleVisitor l;
  l.init(v, diffs);
  t2 = std::chrono::high_resolution_clock::now();
  const double diff_lambda_real_double_setup =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout
    << "LambdaDoubleVisitor setup-time   :"
    << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
    << " us" << std::endl;

  t1 = std::chrono::high_resolution_clock::now();
  for (unsigned j = 0; j < n_runs; j++)
    {
      l.call(res1, vals.data());
      for (unsigned k = 0; k < diffs.size(); k++)
        {
          r2 += res1[k];
        }
    }
  t2 = std::chrono::high_resolution_clock::now();
  const double diff_lambda_real_double_call =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout
    << "LambdaDoubleVisitor run-time     :"
    << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
    << " us" << std::endl;

  deallog << "Value (lambda): " << r2 << std::endl;

#ifdef HAVE_SYMENGINE_LLVM
  // Optimisation via LLVM JIT optimiser

  t1 = std::chrono::high_resolution_clock::now();
  SE::LLVMDoubleVisitor l2;
  l2.init(v, diffs);
  t2 = std::chrono::high_resolution_clock::now();
  std::cout
    << "LLVMDoubleVisitor setup-time     :"
    << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
    << " us" << std::endl;

  t1 = std::chrono::high_resolution_clock::now();
  for (unsigned j = 0; j < 1000; j++)
    {
      l2.call(res2, vals.data());
      for (unsigned k = 0; k < diffs.size(); k++)
        {
          r3 += res2[k];
        }
    }
  t2 = std::chrono::high_resolution_clock::now();
  std::cout
    << "LLVMDoubleVisitor run-time       :"
    << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()
    << " us : Value: " << r3 << std::endl;

  deallog << "Value (llvm): " << r3 << std::endl;
#endif

  deallog << "OK" << std::endl;

  return 0;
}
