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
template BlockVector<double>& BlockVector<double>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<float>(
  const BlockVector<float>&);
template void BlockVector<double>::reinit<double>(const BlockVector<double>&,
						  const bool);
template void BlockVector<double>::reinit<float>(const BlockVector<float>&,
						 const bool);
template void BlockVector<double>::equ<double>(const double,
					       const BlockVector<double>&);
template void BlockVector<double>::equ<float>(const double,
					      const BlockVector<float>&);
template void BlockVector<double>::scale<double>(const BlockVector<double>&);
template void BlockVector<double>::scale<float>(const BlockVector<float>&);

template class BlockVector<float>;
template BlockVector<float>& BlockVector<float>::DEAL_II_MEMBER_OP_TEMPLATE_INST operator=<double>(
  const BlockVector<double>&);
template void BlockVector<float>::reinit<double>(const BlockVector<double>&,
						 const bool);
template void BlockVector<float>::reinit<float>(const BlockVector<float>&,
						const bool);
template void BlockVector<float>::equ<double>(const float,
					      const BlockVector<double>&);
template void BlockVector<float>::equ<float>(const float,
					     const BlockVector<float>&);
template void BlockVector<float>::scale<double>(const BlockVector<double>&);
template void BlockVector<float>::scale<float>(const BlockVector<float>&);

namespace BlockVectorIterators
{
  template class Iterator<double,false>;
  template class Iterator<double,true>;

  template class Iterator<float,false>;
  template class Iterator<float,true>;
}
