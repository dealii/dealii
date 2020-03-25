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

// Check that SymEngine can do delayed substitution of explicitly dependent
// symbolic functions and evaluate their derivatives as well.
// This is inline with what is required to implement some material laws that
// have internal variables

#include "../tests.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra-semi"

#include <symengine/add.h>
#include <symengine/basic.h>
#include <symengine/dict.h>
#include <symengine/eval_double.h>
#include <symengine/functions.h>
#include <symengine/integer.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/real_double.h>
#include <symengine/symbol.h>

#pragma GCC diagnostic pop

#include <iostream>
#include <vector>

namespace SE = SymEngine;

int
main(int argc, char *argv[])
{
  initlog();

  // 0. Define some independent variables
  // Here we have two definitions of Q: One it is treated as a completely
  // indepdent variable, and the other where its dependence on C is
  // explicitly defined
  SE::RCP<const SE::Symbol> C  = SE::symbol("C");  // Primary variable
  SE::RCP<const SE::Symbol> Qi = SE::symbol("Qi"); // Some internal variable
  SE::RCP<const SE::Basic>  Q =
    SE::function_symbol("Q", C); // This sets up the dependence of Q on C

  deallog << *C << ",  " << *Qi << ",  " << *Q << std::endl;

  // 1. Define a function with Q being independent of C
  SE::RCP<const SE::Basic> f_CQ_symb =
    SE::mul(SE::real_double(0.5), SE::mul(Qi, SE::pow(C, SE::integer(2))));

  deallog << "f_CQ_symb: " << *f_CQ_symb << std::endl;

  // 2. Compute the partial derivative of f wrt C. This is a "partial
  // derivative" as Q != Q(C).
  SE::RCP<const SE::Basic> df_CQ_dC_symb = f_CQ_symb->diff(C);

  deallog << "df_CQ_dC_symb: " << *df_CQ_dC_symb << std::endl;

  // 3. Substitute Qi -> Q=Q(C) to now make SymEngine aware of Q's dependence
  // on C.
  SE::map_basic_basic int_var_dict;
  int_var_dict[Qi] = Q;
  SE::RCP<const SE::Basic> df_CQ_dC_symb_subs =
    df_CQ_dC_symb->subs(int_var_dict);

  deallog << "df_CQ_dC_symb_subs: " << *df_CQ_dC_symb_subs << std::endl;

  // 4. Compute total derivative of df_dC wrt C.
  SE::RCP<const SE::Basic> D2f_CQ_DC_dC_symb = df_CQ_dC_symb_subs->diff(C);

  deallog << "D2f_CQ_DC_dC_symb: " << *D2f_CQ_DC_dC_symb << std::endl;

  // 5. Perform some numerical substitutions
  SE::map_basic_basic all_var_dict;
  all_var_dict[C]  = SE::real_double(5);
  all_var_dict[Qi] = SE::real_double(3);
  all_var_dict[Q] =
    SE::mul(SE::real_double(0.6),
            C); // This explicitly defines the relationship between Q and C

  SE::RCP<const SE::Basic> f_CQ_subs = f_CQ_symb->subs(all_var_dict);
  deallog << "f_CQ_subs: " << *f_CQ_subs << std::endl;

  SE::RCP<const SE::Basic> df_CQ_dC_subs = df_CQ_dC_symb->subs(all_var_dict);
  deallog << "df_CQ_dC: subs: " << *df_CQ_dC_subs << std::endl
          << "df_CQ_dC: eval: " << SE::eval_double(*df_CQ_dC_subs) << std::endl;

  // This first substitution should presumably convert Q(C)->0.6C and thus
  // dQ(C)_dC into 0.6
  SE::RCP<const SE::Basic> D2f_CQ_DC_dC_subs_1 =
    D2f_CQ_DC_dC_symb->subs(all_var_dict);
  // This second substitution should fill out the remaining entries for C
  SE::RCP<const SE::Basic> D2f_CQ_DC_dC_subs =
    D2f_CQ_DC_dC_subs_1->subs(all_var_dict);
  deallog << "D2f_CQ_DC_dC: subs 1: " << *D2f_CQ_DC_dC_subs_1 << std::endl
          << "D2f_CQ_DC_dC: subs 2: " << *D2f_CQ_DC_dC_subs << std::endl
          << "D2f_CQ_DC_dC: eval: " << SE::eval_double(*D2f_CQ_DC_dC_subs)
          << std::endl;

  return 0;
}
