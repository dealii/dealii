//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/petsc_block_sparse_matrix.h>


#ifdef DEAL_II_USE_PETSC

namespace PETScWrappers
{
  BlockSparseMatrix::BlockSparseMatrix ()
  {}
  
  

  BlockSparseMatrix::~BlockSparseMatrix ()
  {}



  BlockSparseMatrix &
  BlockSparseMatrix::operator = (const BlockSparseMatrix &m)
  {
    BaseClass::operator = (m);

    return *this;
  }



  BlockSparseMatrix &
  BlockSparseMatrix::operator = (const double d)
  {
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
        this->block(r,c) = d;

    return *this;
  }


  void
  BlockSparseMatrix::
  reinit (const unsigned int n_block_rows,
          const unsigned int n_block_columns)
  {
                                     // first delete previous content of
                                     // the subobjects array
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
        {
          BlockType *p = this->sub_objects[r][c];
          this->sub_objects[r][c] = 0;
          delete p;
        }
  
    this->sub_objects.reinit (0,0);

                                     // then resize. set sizes of blocks to
                                     // zero. user will later have to call
                                     // collect_sizes for this
    this->sub_objects.reinit (n_block_rows,
                              n_block_columns);
    this->row_block_indices.reinit (n_block_rows, 0);
    this->column_block_indices.reinit (n_block_columns, 0);

                                     // and reinitialize the blocks
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      for (unsigned int c=0; c<this->n_block_cols(); ++c)
        {
          BlockType *p = new BlockType();
          this->sub_objects[r][c] = p;
        }
  }



  void
  BlockSparseMatrix::collect_sizes ()
  {
    BaseClass::collect_sizes ();
  }
  
}

  
#endif
