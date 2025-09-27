// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/conditional_ostream.h>

#include <ostream>


DEAL_II_NAMESPACE_OPEN

ConditionalOStream::ConditionalOStream(std::ostream &stream, const bool active)
  : output_stream(stream)
  , active_flag(active)
{}


void
ConditionalOStream::set_condition(bool flag)
{
  active_flag = flag;
}


bool
ConditionalOStream::is_active() const
{
  return active_flag;
}

DEAL_II_NAMESPACE_CLOSE
