//----------------------------  anna_5.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003 by the deal.II authors and Anna Schneebeli
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  anna_5.cc  ---------------------------


// this program used to fail at one point in time.

#include <lac/block_sparsity_pattern.h>
#include <fstream>

int main ()
{
  std::ofstream o("anna_5.output");
  BlockSparsityPattern      sparsity_pattern;
};
