//----------------------------  sparse_matrix.2.templates  ---------------------------
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
//----------------------------  sparse_matrix.2.templates  ---------------------------


// Driver file for SparseMatrix functions with two types.

// TYPEMAT and TYPE2 are defined in sparsematrix?.cc

template SparseMatrix<TYPEMAT> &
SparseMatrix<TYPEMAT>::copy_from<TYPE2> (const SparseMatrix<TYPE2> &);

template 
void SparseMatrix<TYPEMAT>::copy_from<TYPE2> (const FullMatrix<TYPE2> &);

template void SparseMatrix<TYPEMAT>::add_scaled<TYPE2> (const TYPEMAT,
							const SparseMatrix<TYPE2> &);
