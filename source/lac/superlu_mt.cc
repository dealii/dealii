// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/superlu_mt.h>

// add ifdef for superlu-mt
#include "slu_mt_ddefs.h"
//#include "slu_mt_sdefs.h"

#include <typeinfo>

DEAL_II_NAMESPACE_OPEN

SuperLU_MT::SuperLU_MT(const AdditionalData &data)
  : additional_data(data)
{}



SuperLU_MT::~SuperLU_MT()
{
  this->clear();
}



void
SuperLU_MT::clear()
{
	Destroy_SuperMatrix_Store(&this->A);
	Destroy_SuperMatrix_Store(&this->b);
}



void
SuperLU_MT::initialize(const matrix_type &A, const vector_type &b)
{
  Assert(A.m() == A.n(), ExcNotQuadratic());
  Assert(A.m() == b.size(), ExcDimensionMismatch(A.m(), b.size()));
  //   Assert(typeid(number_type) ==
  //   typeid(matrix_type::value_type),ExcMessage("Types mismatch"));

  this->n_rows = A.m();
  this->n_cols = A.n();

  // Actually, the factorize routine of SuperLU-MT only support the
  // compressed column format (SLU-NC). However, the solver routine also
  // accepts the compressed row format (SLU-NR) and then internally converts
  // it into SLU-NC format in order to generate the desired results.
  // Since, dealii stores matrices in a (modified) crs format. We copy
  // the data into the crs format  (SLU-NR) of
  // SuperLU-MT in similar manner to SparseDirectUPFPACK solver.
  this->rowptr.resize(this->n_rows + 1);
  this->colind.resize(A.n_nonzero_elements());
  this->nzval.resize(A.n_nonzero_elements());

  // first fill the rowptr array
  this->rowptr[0] = 0;
  for (size_type row_ctr = 1; row_ctr <= this->n_rows; ++row_ctr)
    this->rowptr[row_ctr] =
      this->rowptr[row_ctr - 1] + A.get_row_length(row_ctr - 1);
  Assert(static_cast<size_type>(this->rowptr.back()) == this->colind.size(),
         ExcInternalError());

  // then copy over matrix elements. note that for sparse matrices,
  // iterators are sorted so that they traverse each row from start to end
  // before moving on to the next row. however, this isn't true for block
  // matrices, so we have to do a bit of book keeping
  {
    // have an array that for each row points to the first entry not yet
    // written to
    std::vector<my_int> row_pointers = this->rowptr;

    // loop over the elements of the matrix row by row, as suggested in the
    // documentation of the sparse matrix iterator class
    for (size_type row = 0; row < A.m(); ++row)
      {
        for (typename matrix_type ::const_iterator p = A.begin(row);
             p != A.end(row);
             ++p)
          {
            // write entry into the first free one for this row
            this->colind[row_pointers[row]] = p->column();
            this->nzval[row_pointers[row]]  = std::real(p->value());

            // then move pointer ahead
            ++row_pointers[row];
          }
      }

    // at the end, we should have written all rows completely
    for (size_type i = 0; i < this->rowptr.size() - 1; ++i)
      Assert(row_pointers[i] == this->rowptr[i + 1], ExcInternalError());
  }

  // now copy the right hand side values
  this->rhs_vals.resize(b.size());
  std::copy(b.begin(), b.end(), this->rhs_vals.begin());

  // now create the \c SuperMatrix objects needed by the solve routine.
  this->create_supermatrices();
}



void
SuperLU_MT::create_supermatrices()
{
  const Stype_t Stype_matrix = Stype_t::SLU_NR;
  const Stype_t Stype_rhs    = Stype_t::SLU_DN;
  const Mtype_t Mtype        = Mtype_t::SLU_GE;

  if (typeid(number_type) == typeid(double))
    {
      const Dtype_t Dtype = Dtype_t::SLU_D;
      dCreate_CompRow_Matrix(&this->A,
                             this->n_rows,
                             this->n_cols,
                             this->nzval.size(),
                             this->nzval.data(),
                             this->colind.data(),
                             this->rowptr.data(),
                             Stype_matrix,
                             Dtype,
                             Mtype);
      dCreate_Dense_Matrix(&this->b,
                           this->n_rows,
                           1,
                           rhs_vals.data(),
                           this->n_rows,
                           Stype_rhs,
                           Dtype,
                           Mtype);
    }
  else if (typeid(number_type) == typeid(float))
    {
      //   const auto& Dtype = Dtype_t::SLU_S;
      //   sCreate_CompCol_Matrix(&A, this->n_rows, this->n_cols,
      //   this->nzval.size(), this->nzval.data(),
      //                          this->colind.data(), this->rowptr.data(),
      //                          Stype_matrix, Dtype, Mtype);
      //   sCreate_Dense_Matrix(&b, this->n_rows, 1, rhs_vals.data(),
      //   this->n_rows, Stype_rhs, Dtype, Mtype);
      Assert(false, ExcInternalError());
    }
  else
    Assert(false, ExcMessage("Type not supported"));
}



void
SuperLU_MT::solve(vector_type &x)
{
  // generate the required datastructures
  int_t              info;
  std::vector<int_t> perm_c(this->n_rows);
  std::vector<int_t> perm_r(this->n_cols);
  SuperMatrix        L;
  SuperMatrix        U;
  get_perm_c(this->additional_data.column_ordering, &this->A, perm_c.data());

  // solve the SLE in parallel
  if (this->A.Dtype == Dtype_t::SLU_D)
    {
      pdgssv(this->additional_data.options_.nprocs,
             &this->A,
             perm_c.data(),
             perm_r.data(),
             &L,
             &U,
             &this->b,
             &info);
    }
  else if (this->A.Dtype == Dtype_t::SLU_S)
    {
      //  psgssv(this->additional_data.options_.nprocs, &this->A, perm_c.data(),
      //  perm_r.data(), &L, &U, &this->b, &info);
      Assert(false, ExcInternalError());
    }

  // copy solution
	DNformat *bStore;
	bStore = (DNformat *)b.Store;
	double *dp;
	dp = (double *)bStore->nzval;
	// \todo check if this manual copy can be replaced by std::copy()
	for (size_type i = 0; i < bStore->lda; ++i)
		x[i] = dp[i];

	// clear the datastructures before leaving
	Destroy_SuperNode_SCP(&L);
	Destroy_CompCol_NCP(&U);
	dp = nullptr;
	bStore = nullptr;
}


DEAL_II_NAMESPACE_CLOSE
