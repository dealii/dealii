//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/block_vector.templates.h>

// explicit instantiations
template class BlockVector<double>;
template class BlockVector<float>;

#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG    
template BlockVector<double>::BlockVector (const BlockVector<float> &);
template BlockVector<float>::BlockVector (const BlockVector<double> &);
#endif

template void BlockVector<double>::reinit<double>(const BlockVector<double>&,
						  const bool);
template void BlockVector<double>::reinit<float>(const BlockVector<float>&,
						 const bool);

template void BlockVector<float>::reinit<double>(const BlockVector<double>&,
						 const bool);
template void BlockVector<float>::reinit<float>(const BlockVector<float>&,
						const bool);

