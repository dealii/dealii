/* $Id$ */

#include <lac/matrix_out.h>


MatrixOut::Options::Options (const bool         show_absolute_values,
			     const unsigned int block_size)
		:
		show_absolute_values (show_absolute_values),
		block_size (block_size)
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
