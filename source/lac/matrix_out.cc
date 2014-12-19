// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#include <deal.II/lac/matrix_out.h>

DEAL_II_NAMESPACE_OPEN


MatrixOut::Options::Options (const bool         show_absolute_values,
                             const unsigned int block_size,
                             const bool         discontinuous)
  :
  show_absolute_values (show_absolute_values),
  block_size (block_size),
  discontinuous (discontinuous)
{}



MatrixOut::~MatrixOut ()
{}



const std::vector<MatrixOut::Patch> &
MatrixOut::get_patches () const
{
  return patches;
}



std::vector<std::string>
MatrixOut::get_dataset_names () const
{
  return std::vector<std::string>(1,name);
}

DEAL_II_NAMESPACE_CLOSE
