//----------------------------  block_vector.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002 by the deal authors
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
template BlockVector<double>& BlockVector<double>::operator=(
  const BlockVector<float>&);
template void BlockVector<double>::reinit(const BlockVector<double>&, bool);
template void BlockVector<double>::reinit(const BlockVector<float>&, bool);
template void BlockVector<double>::equ(double, const BlockVector<double>&);
template void BlockVector<double>::equ(double, const BlockVector<float>&);

template class BlockVector<float>;
template BlockVector<float>& BlockVector<float>::operator=(
  const BlockVector<double>&);
template void BlockVector<float>::reinit(const BlockVector<double>&, bool);
template void BlockVector<float>::reinit(const BlockVector<float>&, bool);
template void BlockVector<float>::equ(float, const BlockVector<double>&);
template void BlockVector<float>::equ(float, const BlockVector<float>&);

namespace BlockVectorIterators
{
  template class Iterator<double,false>;
  template class Iterator<double,true>;

  template class Iterator<float,false>;
  template class Iterator<float,true>;
};
