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
// ORIG-DATE:  2-Jan-03 at 11:02:19 bu Thomas Leurent
//  LAST-MOD:  5-Oct-04 by Jason Kraftcheck
//
// DESCRIPTION:
// ============
/*! \file MsqHessian.hpp

 The MsqHessian class stores a sparse hessian for a given objective 
 function. The objective function must be C2 and such that its hessian
 has non-zero entries only for the duplet of derivatives corresponding 
 to nodes of a same element. 

 \author Thomas Leurent
*/


#ifndef MsqHessian_hpp
#define MsqHessian_hpp

#include "Mesquite.hpp"
#include "Matrix3D.hpp"
#include "PatchData.hpp"
#include "MsqTimer.hpp"

#include <iosfwd>
 
namespace MESQUITE_NS
{
  class ObjectiveFunction;
  
  /*!
    \class MsqHessian
    \brief Vector3D is the object that effeciently stores the objective function
    Hessian each entry is a Matrix3D object (i.e. a vertex Hessian). 
  */
  class MESQUITE_EXPORT MsqHessian
    {
    protected:  // data accessed directly in tests. 
     
      Matrix3D* mEntries;        //!< CSR block entries. size: nb of nonzero blocks, i.e. mRowStart[mSize] . 
      size_t* mRowStart;        //!< start of each row in mEntries. size: nb of vertices (mSize).
      size_t* mColIndex;  //!< CSR block structure: column indexes of the row entries. 

      size_t mSize; //!< number of rows (or number of columns, this is a square matrix).
    
      Matrix3D* mPreconditioner;
      size_t precondArraySize;
    
      Vector3D* mR; //!< array used in the CG solver
      Vector3D* mZ; //!< array used in the CG solver
      Vector3D* mP; //!< array used in the CG solver
      Vector3D* mW; //!< array used in the CG solver
      size_t cgArraySizes; //!< size of arrays allocated in the CG solver.
      size_t maxCGiter; //!< max nb of iterations of the CG solver.
    
    public:
      MsqHessian();
      ~MsqHessian();
    
      void initialize(PatchData &pd, MsqError &err);
      void initialize( const MsqHessian& other );
      
      inline void zero_out();
 
      size_t size() const {return mSize;}
 
      //! returns the diagonal blocks, memory must be allocated before call.
      void get_diagonal_blocks(std::vector<Matrix3D> &diag, MsqError &err) const;

      Matrix3D* get_block(size_t i, size_t j);
      const Matrix3D* get_block(size_t i, size_t j) const;

      //inline void accumulate_entries(PatchData &pd, const size_t &elem_index,
      //                               Matrix3D* const &mat3d_array, MsqError &err);

      void compute_preconditioner(MsqError &err);
      void apply_preconditioner(Vector3D zloc[], Vector3D rloc[], MsqError &err);

      void cg_solver(Vector3D x[], Vector3D b[], MsqError &err);
      //! Hessian - vector product, summed with a second vector (optional).
      friend void axpy(Vector3D res[], size_t size_r,
                       const MsqHessian &H, const Vector3D x[], size_t size_x,
                       const Vector3D y[], size_t size_y, MsqError &err);
      
      //! r = this * x, where r and x are arrays of length size().
      void product( Vector3D* r, const Vector3D* x ) const;

      friend std::ostream& operator<<( std::ostream&, const MsqHessian& );
   
      inline void add( size_t row, size_t col, const Matrix3D& m, MsqError& err );
    
      inline void scale( double value );
      
      void add( const MsqHessian& other ); 
      
        // Release all storage.  Object is invalid until next call
        // to initialize(..).
      void clear();
      
      inline double norm() const; //!< Frobenius norm
      
    private:
      MsqHessian& operator=( const MsqHessian& h );
      MsqHessian( const MsqHessian& copy ); 
    };

  inline double MsqHessian::norm() const {
    double sum = 0.0;
    for (size_t i = 0; i < mRowStart[mSize]; ++i)
      sum += Frobenius_2(mEntries[i]);
    return sqrt(sum);
  }

  /*! Sets all Hessian entries to zero. This is usually used before 
    starting to accumulate elements hessian in the objective function
    hessian. */
  inline void MsqHessian::zero_out()
    {
      if (mSize==0) return; // empty hessian.
    
      size_t i;
      for (i=0; i<mRowStart[mSize]; ++i) {
        mEntries[i].zero();
      }
    }
    
  inline void MsqHessian::scale( double value )
    {
      for (size_t i = 0; i < mRowStart[mSize]; ++i) 
        mEntries[i] *= value;
    }

  
  /*! Accumulates entries of an element hessian into an objective function
    hessian. Make sure to use zero_out() before starting the accumulation
    process. 

    \param pd: PatchData in that contains the element which Hessian
    we are accumulating in the Hessian matrix. This must be the same
    PatchData that was used in MsqHessian::initialize().
    \param elem_index: index of the element in the PatchData.
    \param mat3d_array: This is the upper triangular part of the element Hessian 
    for all nodes, including fixed nodes, for which the entries must be null Matrix3Ds.
    \param nb_mat3d. The size of the mat3d_array: (n+1)n/2, where n is
    the number of nodes in the element.
  */
//  inline void MsqHessian::accumulate_entries(PatchData &pd, const size_t &elem_index,
//                                             Matrix3D* const &mat3d_array, MsqError &err)
//    {
//      if (&pd != origin_pd) {
//        MSQ_SETERR(err)( 
//                    "Cannot accumulate elements from a different patch. "
//                    "Use MsqHessian::initialize first.",
//                    MsqError::INVALID_ARG ); 
//        return;
//      }
//
//      size_t nve = pd.get_element_array(err)[elem_index].vertex_count(); 
//      size_t* ve = pd.get_element_array(err)[elem_index].get_vertex_index_array();
//      size_t e = mAccumElemStart[elem_index];
//      size_t i, j, c = 0;
//      for (i = 0; i < nve; ++i)
//      {
//        for (j = i; j < nve; ++j)
//        {
//          if (ve[i] < mSize && ve[j] < mSize)
//          {
//            int k = mAccumulation[e++];
//            if (k >= 0)
//              mEntries[k] += mat3d_array[c];
//            else
//              mEntries[-k].plus_transpose_equal( mat3d_array[c] );
//          }
//          ++c;
//        }
//      }
//    }
   
  inline void MsqHessian::add( size_t row, size_t col, const Matrix3D& m, MsqError& err )
  {
    if (row <= col) {
      for (size_t i = mRowStart[row]; i != mRowStart[row+1]; ++i)
        if (mColIndex[i] == col) {
          mEntries[i] += m;
          return;
        }
    }
    else {
      for (size_t i = mRowStart[col]; i != mRowStart[col+1]; ++i)
        if (mColIndex[i] == row) {
          mEntries[i].plus_transpose_equal( m );
          return;
        }
    }
    
    MSQ_SETERR(err)(MsqError::INVALID_ARG, 
                   "Hessian entry (%lu,%lu) does not exist.",
                   (unsigned long)row, (unsigned long)col);
  }
   
  /*!
    \param res: array of Vector3D in which the result is stored.
    \param size_r: size of the res array.
    \param x: vector multiplied by the Hessian.
    \param size_x: size of the x array.
    \param y: vector added to the Hessian vector product. Set to 0 (NULL) if not needed.
    \param size_y: size of the y array. Set to 0 if not needed.
  */
  inline void axpy(Vector3D res[], size_t size_r,
                   const MsqHessian &H, const Vector3D x[], size_t size_x,
                   const Vector3D y[], size_t size_y, MsqError &/*err*/)
    {
      if ((size_r != H.mSize) || (size_x != H.mSize) ||
          (size_y != H.mSize && size_y != 0)) {
        // throw an error
      }

      Vector3D tmpx, tmpm; // for cache opt.
      size_t* col = H.mColIndex;
      const size_t nn = H.mSize;
      size_t rl; // row length
      size_t el; // entries index
      size_t lo;
      size_t c;  // column index
      size_t i, j;
     
      if (y != 0) {
        for (i = 0; i < nn; ++i) {
          res[i] = y[i];
        }
      } 
      else {          // y == 0
        for (i = 0; i < nn; ++i) {
          res[i] = 0.;
        }
      }
      
      el = 0;
      for (i = 0; i < nn; ++i) {
        rl = H.mRowStart[i+1] - H.mRowStart[i];
        lo = *col++;

        // Diagonal entry
        tmpx = x[i];
        eqAx(tmpm, H.mEntries[el], tmpx);
        ++el;
        
        //Non-diagonal entries
        for (j = 1; j < rl; ++j) {
          c = *col++;
          //          res[i] += H.mEntries[e] * x[c];
          plusEqAx(tmpm, H.mEntries[el], x[c]);
          //          res[c] += transpose(H.mEntries[e]) * tmpxi;
          plusEqTransAx(res[c], H.mEntries[el], tmpx);
          ++el;
        }
        res[lo] += tmpm;
      }
    }

  /*! Computes \f$ z=M^{-1}r \f$ . */
  inline void MsqHessian::apply_preconditioner(Vector3D zloc[],
                                               Vector3D rloc[],
                                               MsqError& /*err*/)
    {
      size_t m;

      for (m=0; m<mSize; ++m) {
#ifdef DIAGONAL_PRECONDITIONER
        // preconditioner is identity matrix for now.
        zloc[m][0] = mPreconditioner[m][0][0] * rloc[m][0]; 
        zloc[m][1] = mPreconditioner[m][1][1] * rloc[m][1]; 
        zloc[m][2] = mPreconditioner[m][2][2] * rloc[m][2]; 
#else
        // z = inv(L^T) * r
        zloc[m][0] = rloc[m][0];
        zloc[m][1] = rloc[m][1] - mPreconditioner[m][0][1] * zloc[m][0];
        zloc[m][2] = rloc[m][2] - mPreconditioner[m][0][2] * zloc[m][0] - mPreconditioner[m][1][2] * zloc[m][1];

        // z = inv(D) * z
        zloc[m][0] *= mPreconditioner[m][0][0];
        zloc[m][1] *= mPreconditioner[m][1][1];
        zloc[m][2] *= mPreconditioner[m][2][2];

        // z = inv(L) * z
        zloc[m][2] = zloc[m][2];
        zloc[m][1] = zloc[m][1] - mPreconditioner[m][1][2] * zloc[m][2];
        zloc[m][0] = zloc[m][0] - mPreconditioner[m][0][1] * zloc[m][1] - mPreconditioner[m][0][2] * zloc[m][2];
#endif
      }
    }

  
  /*! Returns a pointer to the Matrix3D block at position i,j if it exist. 
    Returns the NULL pointer if position i,j (0-based) is a NULL entry.
    Note that block i,j must be in the upper triangular part of the 
    (symetric) hessian. */
  inline Matrix3D* MsqHessian::get_block(size_t i, size_t j)
    {
      size_t c;
      
      if (i >= mSize || j >= mSize || j < i)
        return NULL;

      for (c=mRowStart[i]; c<mRowStart[i+1]; ++c) {
        if (mColIndex[c] == j)
          return ( mEntries + c );
      }
    
      // if there is no block at position i,j (zero entry).
      return NULL;
    }
  inline const Matrix3D* MsqHessian::get_block(size_t i, size_t j) const
    { return const_cast<MsqHessian*>(this)->get_block(i,j); }

  /* ------------------ I/O ----------------- */

 std::ostream& operator<<(std::ostream &s, const MsqHessian &h);

} // namespace

#endif // MsqHessian_hpp
