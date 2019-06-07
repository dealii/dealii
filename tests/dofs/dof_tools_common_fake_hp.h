// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_dof_tools_common_fake_hp_h
#define dealii_dof_tools_common_fake_hp_h
// fake hp check_this function that does nothing

template <>
void
check_this(const hp::DoFHandler<1, 1> &dof_handler)
{
  // nothing to do here
}

template <>
void
check_this(const hp::DoFHandler<1, 2> &dof_handler)
{
  // nothing to do here
}

template <>
void
check_this(const hp::DoFHandler<1, 3> &dof_handler)
{
  // nothing to do here
}

template <>
void
check_this(const hp::DoFHandler<2, 1> &dof_handler)
{
  // nothing to do here
}

template <>
void
check_this(const hp::DoFHandler<2, 2> &dof_handler)
{
  // nothing to do here
}

template <>
void
check_this(const hp::DoFHandler<2, 3> &dof_handler)
{
  // nothing to do here
}

template <>
void
check_this(const hp::DoFHandler<3, 1> &dof_handler)
{
  // nothing to do here
}

template <>
void
check_this(const hp::DoFHandler<3, 2> &dof_handler)
{
  // nothing to do here
}

template <>
void
check_this(const hp::DoFHandler<3, 3> &dof_handler)
{
  // nothing to do here
}

#endif // dealii_dof_tools_common_fake_hp_h
