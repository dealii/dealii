//----------------------------  mg_base.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_base.cc  ---------------------------


#include <multigrid/mg_base.h>
#include <multigrid/mg_smoother.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <cmath>




 template <class VECTOR>
 MGTransfer<VECTOR>::~MGTransfer()
 {};


 // Explicit instantiations

 template class MGTransfer<Vector<double> >;
 template class MGTransfer<Vector<float> >;
 template class MGTransfer<BlockVector<double> >;
 template class MGTransfer<BlockVector<float> >;
