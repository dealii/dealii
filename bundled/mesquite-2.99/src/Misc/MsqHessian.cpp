/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Todd Munson <tmunson@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tmunson@mcs.anl.gov
//
// ORIG-DATE:  2-Jan-03 at 11:02:19 by Thomas Leurent
//  LAST-MOD: 26-Nov-03 at 15:47:42 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file MsqHessian.cpp

  \author Thomas Leurent

*/

#include "MsqHessian.hpp"
#include "MsqTimer.hpp"


#include <cmath>
#include <iostream>

namespace MESQUITE_NS {

MsqHessian::MsqHessian() :
  mEntries(0), mRowStart(0), mColIndex(0), 
  mSize(0), 
  mPreconditioner(0), precondArraySize(0),
  mR(0), mZ(0), mP(0), mW(0), cgArraySizes(0), maxCGiter(50)
{ }


MsqHessian::~MsqHessian()
{
  clear();
}

void MsqHessian::clear()
{
  mSize = precondArraySize = cgArraySizes = 0;

  delete[] mEntries; mEntries = 0;
  delete[] mRowStart; mRowStart = 0;
  delete[] mColIndex; mColIndex = 0;

  delete[] mPreconditioner; mPreconditioner = 0;

  delete[] mR; mR = 0;
  delete[] mZ; mZ = 0;
  delete[] mP; mP = 0;
  delete[] mW; mW = 0;
}

  
/*! \brief creates a sparse structure for a Hessian, based on the
  connectivity information contained in the PatchData.
  Only the upper triangular part of the Hessian is stored. */
void MsqHessian::initialize(PatchData &pd, MsqError &err)
{
  MSQ_FUNCTION_TIMER( "MsqHession::initialize" );
  delete[] mEntries;
  delete[] mRowStart;
  delete[] mColIndex;
  
  size_t num_vertices = pd.num_free_vertices();
  size_t num_elements = pd.num_elements();
  size_t const * vtx_list;
  size_t e, r, rs, re, c, cs, ce, nz, nnz, nve, i, j;
  MsqMeshEntity* patchElemArray = pd.get_element_array(err); MSQ_CHKERR(err);

  if (num_vertices == 0) {
    MSQ_SETERR( err )( "No vertices in PatchData", MsqError::INVALID_ARG);
    return;
  }

  mSize = num_vertices;

  // Calculate the offsets for a CSC representation of the accumulation
  // pattern.

  size_t* col_start = new size_t[num_vertices + 1];
  //mAccumElemStart = new size_t[num_elements+1];
  //mAccumElemStart[0] = 0;
  
  for (i = 0; i < num_vertices; ++i) {
    col_start[i] = 0;
  }

  for (e = 0; e < num_elements; ++e) {
    nve = patchElemArray[e].node_count();
    vtx_list = patchElemArray[e].get_vertex_index_array();
    int nfe = 0;
    
    for (i = 0; i < nve; ++i) {
      r = vtx_list[i];
      if (r < num_vertices)
        ++nfe;
      
      for (j = i; j < nve; ++j) {
        c = vtx_list[j];

        if (r <= c) {
          if (c < num_vertices)
            ++col_start[c];
        }
        else {
          if (r < num_vertices)
            ++col_start[r];
        }
      }
    }
    //mAccumElemStart[e+1] = mAccumElemStart[e] + (nfe+1)*nfe/2;
  }

  nz = 0;
  for (i = 0; i < num_vertices; ++i) {
    j = col_start[i];
    col_start[i] = nz;
    nz += j;
  }
  col_start[i] = nz;

  // Finished putting matrix into CSC representation

  int* row_instr = new int[5*nz];
  size_t* row_index = new size_t[nz];

  nz = 0;
  for (e = 0; e < num_elements; ++e) {
    nve = patchElemArray[e].node_count();
    vtx_list = patchElemArray[e].get_vertex_index_array();

    for (i = 0; i < nve; ++i) {
      r = vtx_list[i];

      for (j = i; j < nve; ++j) {
        c = vtx_list[j];

        if (r <= c) {
          if (c < num_vertices) {
            row_index[col_start[c]] = r;
            row_instr[col_start[c]] = nz++;
            ++col_start[c];
          }
        }
        else {
          if (r < num_vertices) {
            row_index[col_start[r]] = c;
            //can't use -nz, but can negate row_instr[col_start[r]]
            row_instr[col_start[r]] = nz++;
            row_instr[col_start[r]] = -row_instr[col_start[r]];
            ++col_start[r];
          }
        }
      }
    }
  }

  for (i = num_vertices-1; i > 0; --i) {
    col_start[i+1] = col_start[i];
  }
  col_start[1] = col_start[0];
  col_start[0] = 0;

  //   cout << "col_start: ";
  //   for (int t=0; t<num_vertices+1; ++t)
  //     cout << col_start[t] << " ";
  //   cout << endl;
  //   cout << "row_index: ";
  //   for (int t=0; t<nz; ++t)
  //     cout << row_index[t] << " ";
  //   cout << endl;
  //   cout << "row_instr: ";
  //   for (int t=0; t<nz; ++t)
  //     cout << row_instr[t] << " ";
  //   cout << endl;
  
  
  // Convert CSC to CSR
  // First calculate the offsets in the row

  size_t* row_start = new size_t[num_vertices + 1];
    
  for (i = 0; i < num_vertices; ++i) {
    row_start[i] = 0;
  }

  for (i = 0; i < nz; ++i) {
    ++row_start[row_index[i]];
  }

  nz = 0;
  for (i = 0; i < num_vertices; ++i) {
    j = row_start[i];
    row_start[i] = nz;
    nz += j;
  }
  row_start[i] = nz;
    
  // Now calculate the pattern

  size_t* col_index = new size_t[nz];
  int* col_instr = new int[nz];

  for (i = 0; i < num_vertices; ++i) {
    cs = col_start[i];
    ce = col_start[i+1];

    while(cs < ce) {
      r = row_index[cs];

      col_index[row_start[r]] = i;
      col_instr[row_start[r]] = row_instr[cs];

      ++row_start[r];
      ++cs;
    }
  }

  for (i = num_vertices-1; i > 0; --i) {
    row_start[i+1] = row_start[i];
  }
  row_start[1] = row_start[0];
  row_start[0] = 0;

  delete[] row_index;

  // Now that the matrix is CSR
  // Column indices for each row are sorted

  // Compaction -- count the number of nonzeros
  mRowStart = col_start;   // don't need to reallocate
  //mAccumulation = row_instr;   // don't need to reallocate
  delete [] row_instr;

  for (i = 0; i <= num_vertices; ++i) {
    mRowStart[i] = 0;
  }

  nnz = 0;
  for (i = 0; i < num_vertices; ++i) {
    rs = row_start[i];
    re = row_start[i+1];

    c = num_vertices;
    while(rs < re) {
      if (c != col_index[rs]) {
        // This is an unseen nonzero

        c = col_index[rs];
        ++mRowStart[i];
        ++nnz;
      }

      //if (col_instr[rs] >= 0) {
      //  mAccumulation[col_instr[rs]] = nnz - 1;
      //}
      //else {
      //  mAccumulation[-col_instr[rs]] = 1 - nnz;
      //}
      
      ++rs;
    }
  }

  nnz = 0;
  for (i = 0; i < num_vertices; ++i) {
    j = mRowStart[i];
    mRowStart[i] = nnz;
    nnz += j;
  }
  mRowStart[i] = nnz;

  delete [] col_instr;

  // Fill in the compacted hessian matrix

  mColIndex = new size_t[nnz];

  for (i = 0; i < num_vertices; ++i) {
    rs = row_start[i];
    re = row_start[i+1];

    c = num_vertices;
    while(rs < re) {
      if (c != col_index[rs]) {
        // This is an unseen nonzero

        c = col_index[rs];
        mColIndex[mRowStart[i]] = c;
        mRowStart[i]++;
      }
      ++rs;
    }
  }

  for (i = num_vertices-1; i > 0; --i) {
    mRowStart[i+1] = mRowStart[i];
  }
  mRowStart[1] = mRowStart[0];
  mRowStart[0] = 0;
  
  delete [] row_start;
  delete [] col_index;

  mEntries = new Matrix3D[nnz]; // On Solaris, no initializer allowed for new of an array 
  for (i=0;i<nnz;++i) mEntries[i] = 0.; // so we initialize all entries manually. 

  //origin_pd = &pd;

  return;
}

void MsqHessian::initialize( const MsqHessian& other )
{
  if (!other.mSize) 
  {
    delete[] mEntries;
    delete[] mRowStart;
    delete[] mColIndex;
    mEntries = 0;
    mRowStart = 0;
    mColIndex = 0;
    mSize = 0;
    return;
  }
    
  if (mSize != other.mSize || mRowStart[mSize] != other.mRowStart[mSize])
  {
    delete[] mEntries;
    delete[] mRowStart;
    delete[] mColIndex;
    
    mSize = other.mSize;
    
    mRowStart = new size_t[mSize + 1];
    mEntries = new Matrix3D[other.mRowStart[mSize]];
    mColIndex = new size_t[other.mRowStart[mSize]];
  }
    
  memcpy( mRowStart, other.mRowStart, sizeof(size_t)*(mSize+1) );
  memcpy( mColIndex, other.mColIndex, sizeof(size_t)*mRowStart[mSize] );
}


void MsqHessian::add( const MsqHessian& other )
{
  assert( mSize == other.mSize );
  assert( !memcmp( mRowStart, other.mRowStart, sizeof(size_t)*(mSize+1) ) );
  assert( !memcmp( mColIndex, other.mColIndex, sizeof(size_t)*mRowStart[mSize] ) );
  for (unsigned i = 0; i < mRowStart[mSize]; ++i)
    mEntries[i] += other.mEntries[i];
}
  

/*! \param diag is an STL vector of size MsqHessian::size() . */
void MsqHessian::get_diagonal_blocks(std::vector<Matrix3D> &diag,
                                     MsqError &/*err*/) const
{
  // make sure we have enough memory, so that no reallocation is needed later.
  if (diag.size() != size()) {
    diag.reserve(size());
  }

  for (size_t i=0; i<size(); ++i) {
    diag[i] = mEntries[mRowStart[i]];
  }
}


/*! compute a preconditioner used in the preconditioned conjugate gradient
  algebraic solver. In fact, this computes \f$ M^{-1} \f$ .
*/
void MsqHessian::compute_preconditioner(MsqError &/*err*/)
{
  // reallocates arrays if size of the Hessian has changed too much.
  if (mSize > precondArraySize || mSize < precondArraySize/10 ) {
    delete[] mPreconditioner;
    mPreconditioner = new Matrix3D[mSize];
  }

  Matrix3D* diag_block;
  double sum, tmp;
  size_t m;
  // For each diagonal block, the (inverted) preconditioner is
  // the inverse of the sum of the diagonal entries.
  for (m=0; m<mSize; ++m) {
    diag_block = mEntries + mRowStart[m]; // Gets block at position m,m .

#if !DIAGONAL_PRECONDITIONER
    // calculate LDL^T factorization of the diagonal block
    //  L = [1 pre[0][1] pre[0][2]]
    //      [0 1         pre[1][2]]
    //      [0 0         1        ]
    //  inv(D) = [pre[0][0] 0         0        ]
    //           [0         pre[1][1] 0        ]
    //           [0         0         pre[2][2]]

    //If the first method of calculating a preconditioner fails,
    // use the diagonal method.
    bool use_diag = false;
    
    if (fabs((*diag_block)[0][0]) < DBL_EPSILON) {
      use_diag = true;
    }
    else {
      mPreconditioner[m][0][0] = 1.0 / (*diag_block)[0][0];
      mPreconditioner[m][0][1] = (*diag_block)[0][1] *
          mPreconditioner[m][0][0];
      mPreconditioner[m][0][2] = (*diag_block)[0][2] *
          mPreconditioner[m][0][0];
      sum = ((*diag_block)[1][1] -
             (*diag_block)[0][1] * mPreconditioner[m][0][1]);
      if(fabs(sum)<=DBL_EPSILON)
        use_diag = true;
      else{
        mPreconditioner[m][1][1] = 1.0 / sum;
        
        tmp = (*diag_block)[1][2] - 
            (*diag_block)[0][2] * mPreconditioner[m][0][1];

        mPreconditioner[m][1][2] = mPreconditioner[m][1][1] * tmp;

        sum =  ((*diag_block)[2][2] - 
                (*diag_block)[0][2]*mPreconditioner[m][0][2] - 
                mPreconditioner[m][1][2]*tmp);

        if(fabs(sum)<= DBL_EPSILON)
          use_diag = true;
        else
          mPreconditioner[m][2][2] = 1.0 / sum;
      
      }
      
    }
    if(use_diag)
#endif
    {    
        // Either this is a fixed vertex or the diagonal block is not
        // invertible.  Switch to the diagonal preconditioner in this
        // case.

      sum = (*diag_block)[0][0] + (*diag_block)[1][1] + (*diag_block)[2][2];
      if (fabs(sum) > DBL_EPSILON) 
        sum = 1 / sum;
      else
        sum = 0.0;
      
      mPreconditioner[m][0][0] = sum;
      mPreconditioner[m][0][1] = 0.0;
      mPreconditioner[m][0][2] = 0.0;
      mPreconditioner[m][1][1] = sum;
      mPreconditioner[m][1][2] = 0.0;
      mPreconditioner[m][2][2] = sum;
    }
  }
}


/*! uses the preconditionned conjugate gradient algebraic solver
  to find d in \f$ H * d = -g \f$ .
  \param x : the solution, usually the descent direction d.
  \param b : -b will be the right hand side. Usually b is the gradient.
*/
void MsqHessian::cg_solver(Vector3D x[], Vector3D b[], MsqError &err)
{
  MSQ_FUNCTION_TIMER( "MsqHessian::cg_solver" );
  
  // reallocates arrays if size of the Hessian has changed too much.
  if (mSize > cgArraySizes || mSize < cgArraySizes/10 ) {
    delete[] mR;
    delete[] mZ;
    delete[] mP;
    delete[] mW;
    mR = new Vector3D[mSize];
    mZ = new Vector3D[mSize];
    mP = new Vector3D[mSize];
    mW = new Vector3D[mSize];
    cgArraySizes = mSize;
  }

  size_t i;
  double alpha_, alpha, beta; 
  double cg_tol = 1e-2; // 1e-2 will give a reasonably good solution (~1%). 
  double norm_g = length(b, mSize);
  double norm_r = norm_g;
  double rzm1; // r^T_{k-1} z_{k-1}
  double rzm2; // r^T_{k-2} z_{k-2}

  this->compute_preconditioner(err); MSQ_CHKERR(err); // get M^{-1} for diagonal blocks

  for (i=0; i<mSize; ++i)  x[i] = 0. ;  
  for (i=0; i<mSize; ++i)  mR[i] = -b[i] ;  // r = -b because x_0 = 0 and we solve H*x = -b
  norm_g *= cg_tol;

  this->apply_preconditioner(mZ, mR, err); // solve Mz = r (computes z = M^-1 r)
  for (i=0; i<mSize; ++i)  mP[i] = mZ[i] ; // p_1 = z_0  
  rzm1 = inner(mZ,mR,mSize); // inner product r_{k-1}^T z_{k-1} 
    
  size_t cg_iter = 0;
  while ((norm_r > norm_g) && (cg_iter < maxCGiter)) {
    ++cg_iter;
      
    axpy(mW, mSize, *this, mP, mSize, 0,0,err); // w = A * p_k
      
    alpha_ = inner(mP,mW,mSize); // alpha_ = p_k^T A p_k
    if (alpha_ <= 0.0) {
      if (1 == cg_iter) {
        for (i=0; i<mSize; ++i)  x[i] += mP[i]; // x_{k+1} = x_k + p_{k+1} 
      }
      break; // Newton goes on with this direction of negative curvature 
    }

    alpha = rzm1 / alpha_;
      
    for (i=0; i<mSize; ++i)  x[i] += alpha*mP[i]; // x_{k+1} = x_k + alpha_{k+1} p_{k+1} 
    for (i=0; i<mSize; ++i)  mR[i] -= alpha*mW[i]; // r_{k+1} = r_k - alpha_{k+1} A p_{k+1} 
    norm_r = length(mR, mSize);
 
    this->apply_preconditioner(mZ, mR, err); // solve Mz = r (computes z = M^-1 r)
      
    rzm2 = rzm1;
    rzm1 = inner(mZ,mR,mSize); // inner product r_{k-1}^T z_{k-1} 
    beta = rzm1 / rzm2;
    for (i=0; i<mSize; ++i)  mP[i] = mZ[i] + beta*mP[i]; // p_k = z_{k-1} + Beta_k * p_{k-1}
  }
}

void MsqHessian::product( Vector3D* v, const Vector3D* x ) const
{
    // zero output array
  memset( v, 0, size() * sizeof(*v) );
  
    // for each row
  const Matrix3D* m_iter = mEntries;
  const size_t* c_iter = mColIndex;
  for (size_t r = 0; r < size(); ++r) {
    
      // diagonal entry
    plusEqAx( v[r], *m_iter, x[r] ); 
    ++m_iter; 
    ++c_iter;
    
      // off-diagonal entires
    for (size_t c = mRowStart[r] + 1; c != mRowStart[r+1]; ++c) {
      plusEqAx( v[r], *m_iter, x[*c_iter] );
      plusEqTransAx( v[*c_iter], *m_iter, x[r] );
      ++m_iter;
      ++c_iter;
    }
  }
}

/* ------------------ I/O ----------------- */

//! Prints out the MsqHessian blocks.
std::ostream& operator<<(std::ostream &s, const MsqHessian &h)
{
  size_t i,j;
  s << "MsqHessian of size: " << h.mSize <<"x"<< h.mSize << "\n";
  for (i=0; i<h.mSize; ++i) {
    s << " ROW " << i << " ------------------------\n";
    for (j=h.mRowStart[i]; j<h.mRowStart[i+1]; ++j) {
      s << "   column " << h.mColIndex[j] << " ----\n";
      s << h.mEntries[j]; 
    } 
  }
  return s;
}

} // namespace Mesquite

