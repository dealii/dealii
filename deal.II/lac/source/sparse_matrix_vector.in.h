//---------------------------------------------------------------------------
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
//---------------------------------------------------------------------------


// Driver file for SparseMatrix functions with two types.

// TYPEMAT and TYPEVEC are defined in sparsematrix?.cc

template void SparseMatrix<TYPEMAT>::
vmult<TYPEVEC> (Vector<TYPEVEC> &,
                const Vector<TYPEVEC> &) const;
template void SparseMatrix<TYPEMAT>::
Tvmult<TYPEVEC> (Vector<TYPEVEC> &,
                 const Vector<TYPEVEC> &) const;
template void SparseMatrix<TYPEMAT>::
vmult_add<TYPEVEC> (Vector<TYPEVEC> &,
                    const Vector<TYPEVEC> &) const;
template void SparseMatrix<TYPEMAT>::
Tvmult_add<TYPEVEC> (Vector<TYPEVEC> &,
                     const Vector<TYPEVEC> &) const;

template TYPEVEC
SparseMatrix<TYPEMAT>::
matrix_norm_square<TYPEVEC> (const Vector<TYPEVEC> &) const;

template TYPEVEC
SparseMatrix<TYPEMAT>::
matrix_scalar_product<TYPEVEC> (const Vector<TYPEVEC> &,
                                const Vector<TYPEVEC> &) const;

template TYPEVEC SparseMatrix<TYPEMAT>::
residual<TYPEVEC> (Vector<TYPEVEC> &,
                   const Vector<TYPEVEC> &,
                   const Vector<TYPEVEC> &) const;

template void SparseMatrix<TYPEMAT>::
precondition_SSOR<TYPEVEC> (Vector<TYPEVEC> &,
                            const Vector<TYPEVEC> &,
                            const TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::
precondition_SOR<TYPEVEC> (Vector<TYPEVEC> &,
                           const Vector<TYPEVEC> &,
                           const TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::
precondition_TSOR<TYPEVEC> (Vector<TYPEVEC> &,
                            const Vector<TYPEVEC> &,
                            const TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::
precondition_Jacobi<TYPEVEC> (Vector<TYPEVEC> &,
                              const Vector<TYPEVEC> &,
                              const TYPEMAT) const;

template void SparseMatrix<TYPEMAT>::
SOR<TYPEVEC> (Vector<TYPEVEC> &,
              const TYPEMAT) const;
template void SparseMatrix<TYPEMAT>::
TSOR<TYPEVEC> (Vector<TYPEVEC> &,
               const TYPEMAT) const;
template void SparseMatrix<TYPEMAT>::
SSOR<TYPEVEC> (Vector<TYPEVEC> &,
               const TYPEMAT) const;
template void SparseMatrix<TYPEMAT>::
PSOR<TYPEVEC> (Vector<TYPEVEC> &,
               const std::vector<unsigned int>&,
               const std::vector<unsigned int>&,
               const TYPEMAT) const;
template void SparseMatrix<TYPEMAT>::
TPSOR<TYPEVEC> (Vector<TYPEVEC> &,
                const std::vector<unsigned int>&,
                const std::vector<unsigned int>&,
                const TYPEMAT) const;
template void SparseMatrix<TYPEMAT>::
SOR_step<TYPEVEC> (Vector<TYPEVEC> &,
                   const Vector<TYPEVEC> &,
                   const TYPEMAT) const;
template void SparseMatrix<TYPEMAT>::
TSOR_step<TYPEVEC> (Vector<TYPEVEC> &,
                    const Vector<TYPEVEC> &,
                    const TYPEMAT) const;
template void SparseMatrix<TYPEMAT>::
SSOR_step<TYPEVEC> (Vector<TYPEVEC> &,
                    const Vector<TYPEVEC> &, 
                    const TYPEMAT) const;
