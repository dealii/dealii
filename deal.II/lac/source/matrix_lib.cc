//----------------------------  vector.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vector.cc  ---------------------------


#include <lac/matrix_lib.templates.h>

MeanValueFilter::MeanValueFilter(unsigned int component)
		:
		component(component)
{}

template void MeanValueFilter::filter(Vector<float>&) const;
template void MeanValueFilter::filter(Vector<double>&) const;
template void MeanValueFilter::filter(BlockVector<float>&) const;
template void MeanValueFilter::filter(BlockVector<double>&) const;
template void MeanValueFilter::vmult(Vector<float>&,
				     const Vector<float>&) const;
template void MeanValueFilter::vmult(Vector<double>&,
				     const Vector<double>&) const;
template void MeanValueFilter::vmult(BlockVector<float>&,
				     const BlockVector<float>&) const;
template void MeanValueFilter::vmult(BlockVector<double>&,
				     const BlockVector<double>&) const;

template void MeanValueFilter::vmult_add(Vector<float>&,
					 const Vector<float>&) const;
template void MeanValueFilter::vmult_add(Vector<double>&,
					 const Vector<double>&) const;
template void MeanValueFilter::vmult_add(BlockVector<float>&,
					 const BlockVector<float>&) const;
template void MeanValueFilter::vmult_add(BlockVector<double>&,
					 const BlockVector<double>&) const;

