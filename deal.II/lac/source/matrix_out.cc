/* $Id$ */

#include <lac/matrix_out.h>


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
