//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>

#ifdef DEAL_II_USE_PETSC

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace MPI
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


    void
    BlockSparseMatrix::
    reinit (const unsigned int n_block_rows,
            const unsigned int n_block_columns)
    {
                                       // first delete previous content of
                                       // the subobjects array
      clear ();

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



    const MPI_Comm &
    BlockSparseMatrix::get_mpi_communicator () const
    {
      return block(0,0).get_mpi_communicator();
    }

  }
}

  

DEAL_II_NAMESPACE_CLOSE

#endif
