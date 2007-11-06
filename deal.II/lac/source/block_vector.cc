//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/block_vector.templates.h>

DEAL_II_NAMESPACE_OPEN

// explicit instantiations for real data types
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

// explicit instantiations for complex data types
template class BlockVector<std::complex<double> >;
template class BlockVector<std::complex<float> >;

#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG    
template BlockVector<std::complex<double> >::BlockVector (const BlockVector<std::complex<float> > &);
template BlockVector<std::complex<float> >::BlockVector (const BlockVector<std::complex<double> > &);
#endif

template void BlockVector<std::complex<double> >::reinit<std::complex<double> >(const BlockVector<std::complex<double> >&,
						  const bool);
template void BlockVector<std::complex<double> >::reinit<std::complex<float> >(const BlockVector<std::complex<float> >&,
						 const bool);

template void BlockVector<std::complex<float> >::reinit<std::complex<double> >(const BlockVector<std::complex<double> >&,
						 const bool);
template void BlockVector<std::complex<float> >::reinit<std::complex<float> >(const BlockVector<std::complex<float> >&,
						const bool);


DEAL_II_NAMESPACE_CLOSE
