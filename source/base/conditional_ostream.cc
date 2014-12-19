// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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


#include <deal.II/base/conditional_ostream.h>

DEAL_II_NAMESPACE_OPEN

ConditionalOStream::ConditionalOStream(std::ostream &stream,
                                       const bool    active)
  :
  output_stream (stream),
  active_flag(active)
{}


void ConditionalOStream::set_condition(bool flag)
{
  active_flag = flag;
}


bool ConditionalOStream::is_active() const
{
  return active_flag;
}

DEAL_II_NAMESPACE_CLOSE
