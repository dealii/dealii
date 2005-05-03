//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <numerics/matrices.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>

#ifdef DEAL_II_USE_PETSC
#  include <lac/petsc_parallel_sparse_matrix.h>
#  include <lac/petsc_sparse_matrix.h>
#  include <lac/petsc_parallel_vector.h>
#  include <lac/petsc_vector.h>
#endif

#include <algorithm>


//TODO:[WB] I don't think that the optimized storage of diagonals is needed (GK)
template <typename number>
void
MatrixTools::apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
				    SparseMatrix<number>  &matrix,
				    Vector<number>   &solution,
				    Vector<number>   &right_hand_side,
				    const bool        preserve_symmetry)
{
				   // Require that diagonals are first
				   // in each row
  Assert (matrix.get_sparsity_pattern().optimize_diagonal(),
	  typename SparsityPattern::ExcDiagonalNotOptimized());
  Assert (matrix.n() == right_hand_side.size(),
	  ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
  Assert (matrix.n() == solution.size(),
	  ExcDimensionMismatch(matrix.n(), solution.size()));
				   // if no boundary values are to be applied
				   // simply return
  if (boundary_values.size() == 0)
    return;


  const unsigned int n_dofs = matrix.m();

				   // if a diagonal entry is zero
				   // later, then we use another
				   // number instead. take it to be
				   // the first nonzero diagonal
				   // element of the matrix, or 1 if
				   // there is no such thing
  number first_nonzero_diagonal_entry = 1;
  for (unsigned int i=0; i<n_dofs; ++i)
    if (matrix.diag_element(i) != 0)
      {
	first_nonzero_diagonal_entry = matrix.diag_element(i);
	break;
      };

  
  std::map<unsigned int,double>::const_iterator dof  = boundary_values.begin(),
						endd = boundary_values.end();
  const SparsityPattern    &sparsity    = matrix.get_sparsity_pattern();
  const unsigned int *sparsity_rowstart = sparsity.get_rowstart_indices();
  const unsigned int *sparsity_colnums  = sparsity.get_column_numbers();
  for (; dof != endd; ++dof)
    {
      Assert (dof->first < n_dofs, ExcInternalError());
      
      const unsigned int dof_number = dof->first;
				       // for each boundary dof:
      
				       // set entries of this line
				       // to zero except for the diagonal
				       // entry. Note that the diagonal
				       // entry is always the first one
				       // for square matrices, i.e.
				       // we shall not set
				       // matrix.global_entry(
				       //     sparsity_rowstart[dof.first])
      const unsigned int last = sparsity_rowstart[dof_number+1];
      for (unsigned int j=sparsity_rowstart[dof_number]+1; j<last; ++j)
	matrix.global_entry(j) = 0.;


				       // set right hand side to
				       // wanted value: if main diagonal
				       // entry nonzero, don't touch it
				       // and scale rhs accordingly. If
				       // zero, take the first main
				       // diagonal entry we can find, or
				       // one if no nonzero main diagonal
				       // element exists. Normally, however,
				       // the main diagonal entry should
				       // not be zero.
				       //
				       // store the new rhs entry to make
				       // the gauss step more efficient
      number new_rhs;
      if (matrix.diag_element(dof_number) != 0.0)
        {
          new_rhs = dof->second * matrix.diag_element(dof_number);
          right_hand_side(dof_number) = new_rhs;
        }
      else
	{
	  matrix.set (dof_number, dof_number,
                      first_nonzero_diagonal_entry);
	  new_rhs = dof->second * first_nonzero_diagonal_entry;
          right_hand_side(dof_number) = new_rhs;
	}


				       // if the user wants to have
				       // the symmetry of the matrix
				       // preserved, and if the
				       // sparsity pattern is
				       // symmetric, then do a Gauss
				       // elimination step with the
				       // present row
      if (preserve_symmetry)
	{
					   // store the only nonzero entry
					   // of this line for the Gauss
					   // elimination step
	  const number diagonal_entry = matrix.diag_element(dof_number);
	  
					   // we have to loop over all
					   // rows of the matrix which
					   // have a nonzero entry in
					   // the column which we work
					   // in presently. if the
					   // sparsity pattern is
					   // symmetric, then we can
					   // get the positions of
					   // these rows cheaply by
					   // looking at the nonzero
					   // column numbers of the
					   // present row. we need not
					   // look at the first entry,
					   // since that is the
					   // diagonal element and
					   // thus the present row
	  for (unsigned int j=sparsity_rowstart[dof_number]+1; j<last; ++j)
	    {
	      const unsigned int row = sparsity_colnums[j];

					       // find the position of
					       // element
					       // (row,dof_number)
	      const unsigned int *
		p = std::lower_bound(&sparsity_colnums[sparsity_rowstart[row]+1],
				     &sparsity_colnums[sparsity_rowstart[row+1]],
				     dof_number);

					       // check whether this line has
					       // an entry in the regarding column
					       // (check for ==dof_number and
					       // != next_row, since if
					       // row==dof_number-1, *p is a
					       // past-the-end pointer but points
					       // to dof_number anyway...)
					       //
					       // there should be such an entry!
	      Assert ((*p == dof_number) &&
		      (p != &sparsity_colnums[sparsity_rowstart[row+1]]),
		      ExcInternalError());

	      const unsigned int global_entry
		= (p - &sparsity_colnums[sparsity_rowstart[0]]);
	      
					       // correct right hand side
	      right_hand_side(row) -= matrix.global_entry(global_entry) /
				      diagonal_entry * new_rhs;
	      
					       // set matrix entry to zero
	      matrix.global_entry(global_entry) = 0.;
	    }
	}

				       // preset solution vector
      solution(dof_number) = dof->second;
    }
}



template <typename number>
void
MatrixTools::apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
				    BlockSparseMatrix<number>  &matrix,
				    BlockVector<number>   &solution,
				    BlockVector<number>   &right_hand_side,
				    const bool             preserve_symmetry)
{
  const unsigned int blocks = matrix.n_block_rows();
  
  Assert (matrix.n() == right_hand_side.size(),
	  ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
  Assert (matrix.n() == solution.size(),
	  ExcDimensionMismatch(matrix.n(), solution.size()));
  Assert (matrix.n_block_rows() == matrix.n_block_cols(),
	  ExcNotQuadratic());
  Assert (matrix.get_sparsity_pattern().get_row_indices() == 
	  matrix.get_sparsity_pattern().get_column_indices(),
	  ExcNotQuadratic());
  Assert (matrix.get_sparsity_pattern().get_column_indices() ==
	  solution.get_block_indices (),
	  ExcBlocksDontMatch ());
  Assert (matrix.get_sparsity_pattern().get_row_indices() ==
	  right_hand_side.get_block_indices (),
	  ExcBlocksDontMatch ());

  for (unsigned int i=0; i<blocks; ++i)
    Assert (matrix.block(i,i).get_sparsity_pattern().optimize_diagonal(),
	    SparsityPattern::ExcDiagonalNotOptimized());

  
				   // if no boundary values are to be applied
				   // simply return
  if (boundary_values.size() == 0)
    return;


  const unsigned int n_dofs = matrix.m();

				   // if a diagonal entry is zero
				   // later, then we use another
				   // number instead. take it to be
				   // the first nonzero diagonal
				   // element of the matrix, or 1 if
				   // there is no such thing
  number first_nonzero_diagonal_entry = 0;
  for (unsigned int diag_block=0; diag_block<blocks; ++diag_block)
    {
      for (unsigned int i=0; i<matrix.block(diag_block,diag_block).n(); ++i)
	if (matrix.block(diag_block,diag_block).diag_element(i) != 0)
	  {
	    first_nonzero_diagonal_entry 
	      = matrix.block(diag_block,diag_block).diag_element(i);
	    break;
	  }
				       // check whether we have found
				       // something in the present
				       // block
      if (first_nonzero_diagonal_entry != 0)
	break;
    }
				   // nothing found on all diagonal
				   // blocks? if so, use 1.0 instead
  if (first_nonzero_diagonal_entry == 0)
    first_nonzero_diagonal_entry = 1;
  

  std::map<unsigned int,double>::const_iterator dof  = boundary_values.begin(),
						endd = boundary_values.end();
  const BlockSparsityPattern &
    sparsity_pattern = matrix.get_sparsity_pattern();
  
				   // pointer to the mapping between
				   // global and block indices. since
				   // the row and column mappings are
				   // equal, store a pointer on only
				   // one of them
  const BlockIndices &
    index_mapping = sparsity_pattern.get_column_indices();
  
				   // now loop over all boundary dofs
  for (; dof != endd; ++dof)
    {
      Assert (dof->first < n_dofs, ExcInternalError());

				       // get global index and index
				       // in the block in which this
				       // dof is located
      const unsigned int dof_number = dof->first;
      const std::pair<unsigned int,unsigned int>
	block_index = index_mapping.global_to_local (dof_number);

				       // for each boundary dof:
      
				       // set entries of this line
				       // to zero except for the diagonal
				       // entry. Note that the diagonal
				       // entry is always the first one
				       // for square matrices, i.e.
				       // we shall not set
				       // matrix.global_entry(
				       //     sparsity_rowstart[dof.first])
				       // of the diagonal block
      for (unsigned int block_col=0; block_col<blocks; ++block_col)
	{
	  const SparsityPattern &
	    local_sparsity = sparsity_pattern.block(block_index.first,
						    block_col);

					   // find first and last
					   // entry in the present row
					   // of the present
					   // block. exclude the main
					   // diagonal element, which
					   // is the diagonal element
					   // of a diagonal block,
					   // which must be a square
					   // matrix so the diagonal
					   // element is the first of
					   // this row.
	  const unsigned int 
	    last  = local_sparsity.get_rowstart_indices()[block_index.second+1],
	    first = (block_col == block_index.first ?
		     local_sparsity.get_rowstart_indices()[block_index.second]+1 :
		     local_sparsity.get_rowstart_indices()[block_index.second]);
	  
	  for (unsigned int j=first; j<last; ++j)
	    matrix.block(block_index.first,block_col).global_entry(j) = 0.;
	};
      

				       // set right hand side to
				       // wanted value: if main diagonal
				       // entry nonzero, don't touch it
				       // and scale rhs accordingly. If
				       // zero, take the first main
				       // diagonal entry we can find, or
				       // one if no nonzero main diagonal
				       // element exists. Normally, however,
				       // the main diagonal entry should
				       // not be zero.
				       //
				       // store the new rhs entry to make
				       // the gauss step more efficient
      number new_rhs;
      if (matrix.block(block_index.first, block_index.first)
	  .diag_element(block_index.second) != 0.0)
	new_rhs = dof->second * 
		  matrix.block(block_index.first, block_index.first)
		  .diag_element(block_index.second);
      else
	{
	  matrix.block(block_index.first, block_index.first)
	    .diag_element(block_index.second)
	    = first_nonzero_diagonal_entry;
	  new_rhs = dof->second * first_nonzero_diagonal_entry;
	};
      right_hand_side.block(block_index.first)(block_index.second)
	= new_rhs;


				       // if the user wants to have
				       // the symmetry of the matrix
				       // preserved, and if the
				       // sparsity pattern is
				       // symmetric, then do a Gauss
				       // elimination step with the
				       // present row. this is a
				       // little more complicated for
				       // block matrices.
      if (preserve_symmetry)
	{
					   // store the only nonzero entry
					   // of this line for the Gauss
					   // elimination step
	  const number diagonal_entry 
	    = matrix.block(block_index.first,block_index.first)
	    .diag_element(block_index.second);
	  
					   // we have to loop over all
					   // rows of the matrix which
					   // have a nonzero entry in
					   // the column which we work
					   // in presently. if the
					   // sparsity pattern is
					   // symmetric, then we can
					   // get the positions of
					   // these rows cheaply by
					   // looking at the nonzero
					   // column numbers of the
					   // present row.
					   //
					   // note that if we check
					   // whether row @p{row} in
					   // block (r,c) is non-zero,
					   // then we have to check
					   // for the existence of
					   // column @p{row} in block
					   // (c,r), i.e. of the
					   // transpose block
	  for (unsigned int block_row=0; block_row<blocks; ++block_row)
	    {
					       // get pointers to the
					       // sparsity patterns of
					       // this block and of
					       // the transpose one
	      const SparsityPattern &this_sparsity
		= sparsity_pattern.block (block_row, block_index.first);
	      const SparsityPattern &transpose_sparsity
		= sparsity_pattern.block (block_index.first, block_row);
	      
					       // traverse the row of
					       // the transpose block
					       // to find the
					       // interesting rows in
					       // the present block.
					       // don't use the
					       // diagonal element of
					       // the diagonal block
	      const unsigned int
		first = (block_index.first == block_row ?
			 transpose_sparsity.get_rowstart_indices()[block_index.second]+1 :
			 transpose_sparsity.get_rowstart_indices()[block_index.second]),
		last  = transpose_sparsity.get_rowstart_indices()[block_index.second+1];
	      
	      for (unsigned int j=first; j<last; ++j)
		{
						   // get the number
						   // of the column in
						   // this row in
						   // which a nonzero
						   // entry is. this
						   // is also the row
						   // of the transpose
						   // block which has
						   // an entry in the
						   // interesting row
		  const unsigned int row = transpose_sparsity.get_column_numbers()[j];

						   // find the
						   // position of
						   // element
						   // (row,dof_number)
						   // in this block
						   // (not in the
						   // transpose
						   // one). note that
						   // we have to take
						   // care of special
						   // cases with
						   // square
						   // sub-matrices
		  const unsigned int *p = 0;
		  if (this_sparsity.n_rows() == this_sparsity.n_cols())
		    {
		      if (this_sparsity.get_column_numbers()
			  [this_sparsity.get_rowstart_indices()[row]]
			  ==
			  block_index.second)
			p = &this_sparsity.get_column_numbers()
			    [this_sparsity.get_rowstart_indices()[row]];
		      else
			p = std::lower_bound(&this_sparsity.get_column_numbers()
					     [this_sparsity.get_rowstart_indices()[row]+1],
					     &this_sparsity.get_column_numbers()
					     [this_sparsity.get_rowstart_indices()[row+1]],
					     block_index.second);
		    }
		  else
		    p = std::lower_bound(&this_sparsity.get_column_numbers()
					 [this_sparsity.get_rowstart_indices()[row]],
					 &this_sparsity.get_column_numbers()
					 [this_sparsity.get_rowstart_indices()[row+1]],
					 block_index.second);

						   // check whether this line has
						   // an entry in the regarding column
						   // (check for ==dof_number and
						   // != next_row, since if
						   // row==dof_number-1, *p is a
						   // past-the-end pointer but points
						   // to dof_number anyway...)
						   //
						   // there should be
						   // such an entry!
						   // note, however,
						   // that this
						   // assertion will
						   // fail sometimes
						   // if the sparsity
						   // pattern is not
						   // symmetric!
		  Assert ((*p == block_index.second) &&
			  (p != &this_sparsity.get_column_numbers()
			   [this_sparsity.get_rowstart_indices()[row+1]]),
			  ExcInternalError());
		  
		  const unsigned int global_entry
		    = (p
		       -
		       &this_sparsity.get_column_numbers()
		       [this_sparsity.get_rowstart_indices()[0]]);

						   // correct right hand side
		  right_hand_side.block(block_row)(row)
		    -= matrix.block(block_row,block_index.first).global_entry(global_entry) /
		    diagonal_entry * new_rhs;
		  
						   // set matrix entry to zero
		  matrix.block(block_row,block_index.first).global_entry(global_entry) = 0.;
		}
	    }
	}

				       // preset solution vector
      solution.block(block_index.first)(block_index.second) = dof->second;
    }
}



#ifdef DEAL_II_USE_PETSC

namespace PETScWrappers
{
  template <typename PETScMatrix, typename PETScVector>
  void
  apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
                         PETScMatrix      &matrix,
                         PETScVector      &solution,
                         PETScVector      &right_hand_side,
                         const bool        preserve_symmetry)
  {
    Assert (preserve_symmetry == false, ExcNotImplemented());

    Assert (matrix.n() == right_hand_side.size(),
            ExcDimensionMismatch(matrix.n(), right_hand_side.size()));
    Assert (matrix.n() == solution.size(),
            ExcDimensionMismatch(matrix.n(), solution.size()));

                                     // if no boundary values are to be applied
                                     // simply return
    if (boundary_values.size() == 0)
      return;


                                     // we have to read and write from this
                                     // matrix (in this order). this will only
                                     // work if we compress the matrix first,
                                     // done here. do the same with the other
                                     // objects just to be on the safe side:
    matrix.compress ();
    solution.compress ();
    right_hand_side.compress ();

    const std::pair<unsigned int, unsigned int> local_range
      = matrix.local_range();
    
                                     // determine the first nonzero diagonal
                                     // entry from within the part of the matrix
                                     // that we can see. if we can't find such
                                     // an entry, take one
    PetscScalar average_nonzero_diagonal_entry = 1;
    for (unsigned int i=local_range.first;
         i<local_range.second; ++i)
      if (matrix.diag_element(i) != 0)
        {
          average_nonzero_diagonal_entry = std::fabs(matrix.diag_element(i));
          break;
        }

                                     // iterate over all fixed degrees of
                                     // freedom that are within the local
                                     // range of this matrix. the function is
                                     // pretty awkward to implement, since we
                                     // can't freely mix reading and writing
                                     // from the matrix without global
                                     // synchronisation, so we first determine
                                     // which are the entries we have to work
                                     // on in a first step, and then come back
                                     // and do that all at once
    {
      matrix.compress ();
      
      std::map<unsigned int,double>::const_iterator
        dof  = boundary_values.begin(),
        endd = boundary_values.end();
      std::vector<std::pair<unsigned int, unsigned int> >
        set_to_zero_entries;
      
      for (; dof != endd; ++dof)
        if ((dof->first >= local_range.first) &&
            (dof->first < local_range.second))
          {
            const unsigned int dof_number = dof->first;
      
                                             // for each constrained dof:
      
                                             // store which entries of this line
                                             // to set to zero except for the
                                             // diagonal entry.
            PETScWrappers::MatrixBase::const_iterator
              p = matrix.begin(dof_number),
              e = matrix.end(dof_number);

            for (; p!=e; ++p)
              {
                Assert (p->row() == dof_number, ExcInternalError());
                
                if (p->column() != dof_number)
                  set_to_zero_entries.push_back (std::make_pair (dof_number,
                                                                 p->column()));
              }
          }

                                       // now set all these entries to zero in
                                       // one bulk operation that requires
                                       // only a single synchronisation:
      matrix.compress ();
      for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
             i = set_to_zero_entries.begin();
           i != set_to_zero_entries.end(); ++i)
        matrix.set (i->first, i->second, 0.);
    }
      

                                           // the next thing is to set right
                                           // hand side to wanted value. note
                                           // that for petsc matrices
                                           // interleaving read with write
                                           // operations is very
                                           // expensive. thus, we here always
                                           // replace the diagonal element,
                                           // rather than first checking
                                           // whether it is nonzero and in
                                           // that case preserving it. this is
                                           // different from the case of
                                           // deal.II sparse matrices treated
                                           // in the other functions.
    {
      matrix.compress ();

      std::map<unsigned int,double>::const_iterator
        dof  = boundary_values.begin(),
        endd = boundary_values.end();

      for (; dof != endd; ++dof)
        if ((dof->first >= local_range.first) &&
            (dof->first < local_range.second))
          {  
            const unsigned int dof_number = dof->first;

            matrix.set (dof_number, dof_number,
                        average_nonzero_diagonal_entry);
            right_hand_side(dof_number)
              = dof->second * average_nonzero_diagonal_entry;

                                             // preset solution vector
            solution(dof_number) = dof->second;
          }
    }
    
    matrix.compress ();
    solution.compress ();
    right_hand_side.compress ();
  }
}



void
MatrixTools::
apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
                       PETScWrappers::SparseMatrix   &matrix,
                       PETScWrappers::Vector   &solution,
                       PETScWrappers::Vector   &right_hand_side,
                       const bool        preserve_symmetry)
{
                                   // simply redirect to the generic function
                                   // used for both petsc matrix types
  PETScWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                        right_hand_side, preserve_symmetry);
}



void
MatrixTools::
apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
                       PETScWrappers::MPI::SparseMatrix   &matrix,
                       PETScWrappers::MPI::Vector   &solution,
                       PETScWrappers::MPI::Vector   &right_hand_side,
                       const bool        preserve_symmetry)
{
                                   // simply redirect to the generic function
                                   // used for both petsc matrix types
  PETScWrappers::apply_boundary_values (boundary_values, matrix, solution,
                                        right_hand_side, preserve_symmetry);
}


#endif


void
MatrixTools::
local_apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
                             const std::vector<unsigned int> &local_dof_indices,
                             FullMatrix<double> &local_matrix,
                             Vector<double>     &local_rhs,
                             const bool          eliminate_columns)
{
  Assert (local_dof_indices.size() == local_matrix.m(),
          ExcDimensionMismatch(local_dof_indices.size(),
                               local_matrix.m()));
  Assert (local_dof_indices.size() == local_matrix.n(),
          ExcDimensionMismatch(local_dof_indices.size(),
                               local_matrix.n()));
  Assert (local_dof_indices.size() == local_rhs.size(),
          ExcDimensionMismatch(local_dof_indices.size(),
                               local_rhs.size()));

                                   // if there is nothing to do, then exit
                                   // right away
  if (boundary_values.size() == 0)
    return;

                                   // otherwise traverse all the dofs used in
                                   // the local matrices and vectors and see
                                   // what's there to do

                                   // if we need to treat an entry, then we
                                   // set the diagonal entry to its absolute
                                   // value. if it is zero, we used to set it
                                   // to one, which is a really terrible
                                   // choice that can lead to hours of
                                   // searching for bugs in programs (I
                                   // experienced this :-( ) if the matrix
                                   // entries are otherwise very large. this
                                   // is so since iterative solvers would
                                   // simply not correct boundary nodes for
                                   // their correct values since the residual
                                   // contributions of their rows of the
                                   // linear system is almost zero if the
                                   // diagonal entry is one. thus, set it to
                                   // the average absolute value of the
                                   // nonzero diagonal elements.
                                   //
                                   // we only compute this value lazily the
                                   // first time we need it.
  double average_diagonal = 0;
  const unsigned int n_local_dofs = local_dof_indices.size();
  for (unsigned int i=0; i<n_local_dofs; ++i)
    {
      const std::map<unsigned int, double>::const_iterator
        boundary_value = boundary_values.find (local_dof_indices[i]);
      if (boundary_value != boundary_values.end())
        {
                                           // remove this row, except for the
                                           // diagonal element
          for (unsigned j=0; j<n_local_dofs; ++j)
            if (i != j)
              local_matrix(i,j) = 0;

                                           // replace diagonal entry by its
                                           // absolute value to make sure that
                                           // everything remains positive, or
                                           // by the average diagonal value if
                                           // zero
          if (local_matrix(i,i) == 0.)
            {
                                               // if average diagonal hasn't
                                               // yet been computed, do so now
              if (average_diagonal == 0.)
                {
                  unsigned int nonzero_diagonals = 0;
                  for (unsigned int k=0; k<n_local_dofs; ++k)
                    if (local_matrix(k,k) != 0.)
                      {
                        average_diagonal += std::fabs(local_matrix(k,k));
                        ++nonzero_diagonals;
                      }
                  if (nonzero_diagonals != 0)
                    average_diagonal /= nonzero_diagonals;
                  else
                    average_diagonal = 0;
                }

                                               // only if all diagonal entries
                                               // are zero, then resort to the
                                               // last measure: choose one
              if (average_diagonal == 0.)
                average_diagonal = 1.;

              local_matrix(i,i) = average_diagonal;
            }
          else
            local_matrix(i,i) = std::fabs(local_matrix(i,i));
          
                                           // and replace rhs entry by correct
                                           // value
          local_rhs(i) = local_matrix(i,i) * boundary_value->second;

                                           // finally do the elimination step
                                           // if requested
          if (eliminate_columns == true)
            {
              for (unsigned int row=0; row<n_local_dofs; ++row)
                if (row != i)
                  {
                    local_rhs(row) -= local_matrix(row,i) * boundary_value->second;
                    local_matrix(row,i) = 0;
                  }
            }
        }
    }
}





// explicit instantiations

template
void
MatrixTools::apply_boundary_values<double> (const std::map<unsigned int,double> &boundary_values,
					    SparseMatrix<double>  &matrix,
					    Vector<double>   &solution,
					    Vector<double>   &right_hand_side,
					    const bool        preserve_symmetry);
template
void
MatrixTools::apply_boundary_values<float> (const std::map<unsigned int,double> &boundary_values,
					   SparseMatrix<float>  &matrix,
					   Vector<float>   &solution,
					   Vector<float>   &right_hand_side,
					   const bool        preserve_symmetry);

template
void
MatrixTools::apply_boundary_values<double> (const std::map<unsigned int,double> &boundary_values,
					    BlockSparseMatrix<double>  &matrix,
					    BlockVector<double>   &solution,
					    BlockVector<double>   &right_hand_side,
					    const bool        preserve_symmetry);
template
void
MatrixTools::apply_boundary_values<float> (const std::map<unsigned int,double> &boundary_values,
					   BlockSparseMatrix<float>  &matrix,
					   BlockVector<float>   &solution,
					   BlockVector<float>   &right_hand_side,
					   const bool        preserve_symmetry);
