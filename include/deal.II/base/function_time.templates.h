// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 1999 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_function_time_templates_h
#define dealii_function_time_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/function_time.h>

DEAL_II_NAMESPACE_OPEN


template <typename Number>
FunctionTime<Number>::FunctionTime(const Number initial_time)
  : time(initial_time)
{}



template <typename Number>
void
FunctionTime<Number>::set_time(const Number new_time)
{
  time = new_time;
}


template <typename Number>
void
FunctionTime<Number>::advance_time(const Number delta_t)
{
  set_time(time + delta_t);
}


DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_function_time_templates_h */
