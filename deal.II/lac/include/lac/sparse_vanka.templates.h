//----------------------------  sparse_vanka.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_vanka.templates.h  ---------------------------
#ifndef __deal2__sparse_vanka_templates_h
#define __deal2__sparse_vanka_templates_h


#include <lac/sparse_vanka.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>

#ifdef DEAL_II_USE_MT
#  include <base/thread_management.h>
#endif

#include <algorithm>
#include <map>


template<typename number>
SparseVanka<number>::SparseVanka(const SparseMatrix<number> &M,
				 const vector<bool>         &selected,
				 const bool                  conserve_mem,
				 const unsigned int          n_threads)
		:
		matrix (&M),
		conserve_mem (conserve_mem),
		selected (selected),
		n_threads (n_threads),
		inverses (M.m(), 0)
{
  Assert (M.m() == M.n(), ExcMatrixNotSquare ());
  Assert (M.m() == selected.size(), ExcInvalidVectorSize(M.m(), selected.size()));

  if (conserve_mem == false)
    compute_inverses ();
}


template<typename number>
SparseVanka<number>::~SparseVanka()
{
  vector<SmartPointer<FullMatrix<float> > >::iterator i;
  for(i=inverses.begin(); i!=inverses.end(); ++i)
    {
      FullMatrix<float> *p = *i;
      *i = 0;
      if (p != 0) delete p;
    }
}


template <typename number>
void
SparseVanka<number>::compute_inverses () 
{
#ifdef DEAL_II_USE_MT
  const unsigned int n_inverses = count (selected.begin(),
					 selected.end(),
					 true);

  const unsigned int n_inverses_per_thread = max(n_inverses / n_threads,
						 1U);
  
				   // set up start and end index for
				   // each of the threads. note that
				   // we have to work somewhat to get
				   // this appropriate, since the
				   // indices for which inverses have
				   // to be computed may not be evenly
				   // distributed in the vector. as an
				   // extreme example consider
				   // numbering of DoFs by component,
				   // then all indices for which we
				   // have to do work will be
				   // consecutive, with other
				   // consecutive regions where we do
				   // not have to do something
  vector<pair<unsigned int, unsigned int> > blocking (n_threads);

  unsigned int c       = 0;
  unsigned int thread  = 0;
  blocking[0].first = 0;
  
  for (unsigned int i=0; (i<matrix->m()) && (thread+1<n_threads); ++i)
    {
      if (selected[i] == true)
	++c;
      if (c == n_inverses_per_thread)
	{
	  blocking[thread].second  = i;
	  blocking[thread+1].first = i;
	  ++thread;

	  c = 0;
	};
    };
  blocking[n_threads-1].second = matrix->m();

  typedef void (SparseVanka<number>::*FunPtr)(unsigned int, unsigned int);
  FunPtr fun_ptr = &SparseVanka<number>::compute_inverses;
  
				   // Now spawn the threads
  Threads::ThreadManager thread_manager;
  for (unsigned int i=0; i<n_threads; ++i)
    Threads::spawn (thread_manager,
		    Threads::encapsulate (fun_ptr)
		    .collect_args (this, blocking[i].first, blocking[i].second));
  
  thread_manager.wait ();
  
#else
  compute_inverses (0, matrix->m());
#endif
};


template <typename number>
void
SparseVanka<number>::compute_inverses (const unsigned int begin,
				       const unsigned int end) 
{
				   // set-up the vector that will be used
				   // by the functions which we call
				   // below.
  vector<unsigned int> local_indices;

				   // traverse all rows of the matrix
				   // which are selected
  for (unsigned int row=begin; row<end; ++row)
    if (selected[row] == true)
      compute_inverse (row, local_indices);
};



template <typename number>
void
SparseVanka<number>::compute_inverse (const unsigned int    row,
				      vector<unsigned int> &local_indices) 
{
				   // first define an alias to the sparsity
				   // pattern of the matrix, since this
				   // will be used quite often
  const SparsityPattern &structure
    = matrix->get_sparsity_pattern();

  const unsigned int row_length = structure.row_length(row);
  inverses[row] = new FullMatrix<float> (row_length,
					 row_length); 

				   // collect the dofs that couple
				   // with @p{row}
  local_indices.resize (row_length);
  for (unsigned int i=0; i<row_length; ++i)
    local_indices[i] = structure.column_number(row, i);
  
				   // Build local matrix
  for (unsigned int i=0; i<row_length; ++i)
    for (unsigned int j=0; j<row_length; ++j)
      {
					 // if DoFs local_index[i] and
					 // local_index[j] couple with
					 // each other, then get the
					 // value from the global
					 // matrix. if not, then leave
					 // the value in the small
					 // matrix at zero
	const unsigned int global_entry =
	  structure(local_indices[i], local_indices[j]);
	if (global_entry != SparsityPattern::invalid_entry)
	  (*inverses[row])(i,j) = matrix->global_entry(global_entry);
      };
  
				   // Compute inverse
  inverses[row]->gauss_jordan();
};

//TODO: implement Tvmult

template<typename number>
template<typename number2>
void
SparseVanka<number>::vmult (Vector<number2>       &dst,
			    const Vector<number2> &src) const
{
				   // first set output vector to zero
  dst.clear ();
				   // then pass on to the function
				   // that actually does the work
  apply_preconditioner (dst, src);
};



template<typename number>
template<typename number2>
void
SparseVanka<number>::apply_preconditioner (Vector<number2>       &dst,
					   const Vector<number2> &src,
					   const vector<bool>    *const dof_mask) const
{
  Assert (dst.size() == src.size(),
	  ExcInvalidVectorSize(dst.size(), src.size()));
  Assert (dst.size() == matrix->m(),
	  ExcInvalidVectorSize(dst.size(), src.size()));
  
				   // first define an alias to the sparsity
				   // pattern of the matrix, since this
				   // will be used quite often
  const SparsityPattern &structure
    = matrix->get_sparsity_pattern();


				   // store whether we shall work on
				   // the whole matrix, or only on
				   // blocks. this variable is used to
				   // optimize access to vectors a
				   // little bit.
  const bool range_is_restricted = (dof_mask != 0);
  
				   // space to be used for local
				   // systems. allocate as much memory
				   // as is the maximum. this
				   // eliminates the need to
				   // re-allocate memory inside the
				   // loop.
  FullMatrix<float> local_matrix (structure.max_entries_per_row(),
				  structure.max_entries_per_row());
  Vector<float> b (structure.max_entries_per_row());
  Vector<float> x (structure.max_entries_per_row());
  
  map<unsigned int, unsigned int> local_index;

				   // traverse all rows of the matrix
				   // which are selected
  const unsigned int n = matrix->m();
  for (unsigned int row=0; row<n; ++row)
    if ((selected[row] == true) &&
	((range_is_restricted == false) || ((*dof_mask)[row] == true)))
      {
	const unsigned int row_length = structure.row_length(row);
	
					 // if we don't store the
					 // inverse matrices, then alias
					 // the entry in the global
					 // vector to the local matrix
					 // to be used
	if (conserve_mem == true)
	  {
	    inverses[row] = &local_matrix;
	    inverses[row]->reinit (row_length, row_length);
	  };
	
	b.reinit (row_length);
	x.reinit (row_length);
					 // mapping between:
					 // 1 column number of all
					 //   entries in this row, and
					 // 2 the position within this
					 //   row (as stored in the
					 //   SparsityPattern object
					 //
					 // since we do not explicitely
					 // consider nonsysmmetric sparsity
					 // patterns, the first element
					 // of each entry simply denotes
					 // all degrees of freedom that
					 // couple with @p{row}.
	local_index.clear ();
	for (unsigned int i=0; i<row_length; ++i)
	  local_index.insert(pair<unsigned int, unsigned int>
			     (structure.column_number(row, i), i));
	
					 // Build local matrix and rhs
	for (map<unsigned int, unsigned int>::const_iterator is=local_index.begin();
	     is!=local_index.end(); ++is)
	  {
					     // irow loops over all DoFs that
					     // couple with the present DoF
	    const unsigned int irow = is->first;
					     // index of DoF irow in the matrix
					     // row corresponding to DoF @p{row}.
					     // runs between 0 and row_length
	    const unsigned int i = is->second;
					     // number of DoFs coupling to
					     // irow (including irow itself)
	    const unsigned int irow_length = structure.row_length(irow);
	    
					     // copy rhs
	    b(i) = src(irow);
	    
					     // for all the DoFs that irow
					     // couples with
	    for (unsigned int j=0; j<irow_length; ++j)
	      {
						 // col is the number of
						 // this dof
		const unsigned int col = structure.column_number(irow, j);
						 // find out whether this DoF
						 // (that couples with @p{irow},
						 // which itself couples with
						 // @p{row}) also couples with
						 // @p{row}.
		const map<unsigned int, unsigned int>::const_iterator js
		  = local_index.find(col);
						 // if not, then still use
						 // this dof to modify the rhs
						 //
						 // note that if so, we already
						 // have copied the entry above
		if (js == local_index.end())
		  {
		    if (!range_is_restricted ||
			((*dof_mask)[col] == true))
		      b(i) -= matrix->raw_entry(irow,j) * dst(col);
		  }
		else
						     // if so, then build the
						     // matrix out of it
		  if (conserve_mem == true)
		    (*inverses[row])(i,js->second) = matrix->raw_entry(irow,j);
	      };
	  };
	
					 // Compute new values
	if (conserve_mem == true)
	  inverses[row]->gauss_jordan();

					 // apply preconditioner
	inverses[row]->vmult(x,b);
	
					 // Distribute new values
	for (map<unsigned int, unsigned int>::const_iterator is=local_index.begin();
	     is!=local_index.end(); ++is)
	  {
	    const unsigned int irow = is->first;
	    const unsigned int i = is->second;

	    if (!range_is_restricted ||
		((*dof_mask)[irow] == true))
	      dst(irow) = x(i);
					       // do nothing if not in
					       // the range
	  };
	
					 // if we don't store the
					 // inverses, then unalias the
					 // local matrix
	if (conserve_mem == true)
	  inverses[row] = 0;
      };
};


template <typename number>
SparseBlockVanka<number>::SparseBlockVanka (const SparseMatrix<number> &M,
					    const vector<bool>         &selected,
					    const unsigned int          n_blocks,
					    const BlockingStrategy      blocking_strategy,
					    const bool                  conserve_memory,
					    const unsigned int          n_threads)
		:
		SparseVanka<number> (M, selected, conserve_memory, n_threads),
                n_blocks (n_blocks),
                dof_masks (n_blocks,
		           vector<bool>(M.m(), false))
{
  compute_dof_masks (M, selected, blocking_strategy);
};


template <typename number>
void
SparseBlockVanka<number>::compute_dof_masks (const SparseMatrix<number> &M,
					     const vector<bool>         &selected,
					     const BlockingStrategy      blocking_strategy) 
{
  Assert (n_blocks > 0, ExcInternalError());

  const unsigned int n_inverses = count (selected.begin(),
					 selected.end(),
					 true);
  
  const unsigned int n_inverses_per_block = max(n_inverses / n_blocks,
						1U);
  
				   // precompute the splitting points
  vector<pair<unsigned int, unsigned int> > intervals (n_blocks);
    
				   // set up start and end index for
				   // each of the blocks. note that
				   // we have to work somewhat to get
				   // this appropriate, since the
				   // indices for which inverses have
				   // to be computed may not be evenly
				   // distributed in the vector. as an
				   // extreme example consider
				   // numbering of DoFs by component,
				   // then all indices for which we
				   // have to do work will be
				   // consecutive, with other
				   // consecutive regions where we do
				   // not have to do something
  if (true)
    {
      unsigned int c       = 0;
      unsigned int block   = 0;
      intervals[0].first   = 0;
      
      for (unsigned int i=0; (i<M.m()) && (block+1<n_blocks); ++i)
	{
	  if (selected[i] == true)
	    ++c;
	  if (c == n_inverses_per_block)
	    {
	      intervals[block].second  = i;
	      intervals[block+1].first = i;
	      ++block;
	      
	      c = 0;
	    };
	};
      intervals[n_blocks-1].second = M.m();
    };

				   // now transfer the knowledge on
				   // the splitting points into the
				   // vector<bool>s that the base
				   // class wants to see. the way how
				   // we do this depends on the
				   // requested blocking strategy
  switch (blocking_strategy)
    {
      case index_intervals:
      {
	for (unsigned int block=0; block<n_blocks; ++block)
	  fill_n (dof_masks[block].begin()+intervals[block].first,
		  intervals[block].second - intervals[block].first,
		  true);
	break;
      };

      case adaptive:
      {
					 // the splitting points for
					 // the DoF have been computed
					 // above already, but we will
					 // only use them to split the
					 // Lagrange dofs into
					 // blocks. splitting the
					 // remaining dofs will be
					 // done now.
	
					 // first count how often the
					 // Lagrange dofs of each
					 // block access the different
					 // dofs
	vector<vector<unsigned int> > access_count (n_blocks,
						    vector<unsigned int>(M.m(), 0));
	
					 // set-up the map that will
					 // be used to store the
					 // indices each Lagrange dof
					 // accesses
	map<unsigned int, unsigned int> local_index;
	const SparsityPattern &structure = M.get_sparsity_pattern();

	for (unsigned int row=0; row<M.m(); ++row)
	  if (selected[row] == true)
	    {
					       // first find out to
					       // which block the
					       // present row belongs
	      unsigned int block_number = 0;
	      while (row>=intervals[block_number].second)
		++block_number;
	      Assert (block_number < n_blocks, ExcInternalError());
	      
					       // now traverse the
					       // matrix structure to
					       // find out to which
					       // dofs number the
					       // present index wants
					       // to write
	      const unsigned int row_length = structure.row_length(row);
	      for (unsigned int i=0; i<row_length; ++i)
		++access_count[block_number][structure.column_number(row, i)];
	    };

					 // now we know that block @p{i}
					 // wants to write to DoF @p{j}
					 // as often as
					 // @p{access_count[i][j]}
					 // times. Let @p{j} be allotted
					 // to the block which
					 // accesses it most often.
					 //
					 // if it is a Lagrange dof,
					 // the of course we leave it
					 // to the block we put it
					 // into in the first place
	for (unsigned int row=0; row<M.m(); ++row)
	  if (selected[row] == true)
	    {
	      unsigned int block_number = 0;
	      while (row>=intervals[block_number].second)
		++block_number;

	      dof_masks[block_number][row] = true;
	    }
	  else
	    {
					       // find out which block
					       // accesses this dof
					       // the most often
	      unsigned int max_accesses     = 0;
	      unsigned int max_access_block = 0;
	      for (unsigned int block=0; block<n_blocks; ++block)
		if (access_count[block][row] > max_accesses)
		  {
		    max_accesses = access_count[block][row];
		    max_access_block = block;
		  };
	      dof_masks[max_access_block][row] = true;
	    };

	break;
      };
       
      default:
	    Assert (false, ExcInternalError());
    };
};

//TODO: implement Tvmult

template <typename number>
template <typename number2>
void SparseBlockVanka<number>::vmult (Vector<number2>       &dst,
				      const Vector<number2> &src) const
{
  dst.clear ();

				   // if no blocking is required, pass
				   // down to the underlying class
  if (n_blocks == 1)
    apply_preconditioner (dst, src);
  else
				     // otherwise: blocking requested
    {
#ifdef DEAL_II_USE_MT
      Threads::ThreadManager thread_manager;
      for (unsigned int block=0; block<n_blocks; ++block)
	Threads::spawn (thread_manager,
			Threads::encapsulate (&SparseVanka<number>::
					      template apply_preconditioner<number2>)
			.collect_args (this, dst, src, &dof_masks[block]));

      thread_manager.wait ();
#else
      for (unsigned int block=0; block<n_blocks; ++block)
	apply_preconditioner (dst, src,
			      &dof_masks[block]);
#endif
    }
}

#endif
