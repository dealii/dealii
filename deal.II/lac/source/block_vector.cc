//----------------------------  block_vector.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_vector.cc  ---------------------------

#include <lac/block_vector.templates.h>

// explicit instantiations
template class BlockVector<double>;
template class BlockVector<float>;

namespace BlockVectorIterators
{
  template class Iterator<double,false>;
  template class Iterator<double,true>;

  template class Iterator<float,false>;
  template class Iterator<float,true>;
};

