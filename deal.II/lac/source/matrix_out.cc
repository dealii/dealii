//----------------------------  matrix_out.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  matrix_out.cc  ---------------------------


#include <lac/matrix_out.h>


MatrixOut::Options::Options (const bool         show_absolute_values,
			     const unsigned int block_size,
			     const bool         discontinuous)
		:
		show_absolute_values (show_absolute_values),
		block_size (block_size),
		discontinuous (discontinuous)
{};



MatrixOut::~MatrixOut () 
{};



const std::vector<MatrixOut::Patch> &
MatrixOut::get_patches () const
{
  return patches;
};



std::vector<std::string> 
MatrixOut::get_dataset_names () const
{
  return std::vector<std::string>(1,name);
};
