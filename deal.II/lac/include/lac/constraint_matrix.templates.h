//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_constraints_templates_h
#define __deal2__dof_constraints_templates_h


#include <base/config.h>
#include <lac/constraint_matrix.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_block_sparse_matrix.h>
#include <lac/petsc_block_sparse_matrix.h>
#include <lac/petsc_parallel_block_sparse_matrix.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/block_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN


template<typename number>
void
ConstraintMatrix::condense (const SparseMatrix<number> &uncondensed,
			    SparseMatrix<number>       &condensed) const
{
				   // create two dummy vectors and enter the
				   // other function
  Vector<number> dummy (0);
  condense (uncondensed, dummy, condensed, dummy);
}



template<typename number>
void
ConstraintMatrix::condense (SparseMatrix<number> &uncondensed) const
{
  Vector<number> dummy (0);
  condense (uncondensed, dummy);
}



template <typename number>
void
ConstraintMatrix::condense (BlockSparseMatrix<number> &uncondensed) const
{
  BlockVector<number> dummy (0);
  condense (uncondensed, dummy);
}



template<class VectorType>
void
ConstraintMatrix::condense (const VectorType &uncondensed,
			    VectorType       &condensed) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (condensed.size()+n_constraints() == uncondensed.size(),
	  ExcDimensionMismatch(condensed.size()+n_constraints(),
			       uncondensed.size()));
  
				   // store for each line of the
				   // vector its new line number after
				   // compression. If the shift is -1,
				   // this line will be condensed away
  std::vector<int> new_line;

  new_line.reserve (uncondensed.size());

  std::vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                                shift           = 0;
  unsigned int n_rows = uncondensed.size();

  if (next_constraint == lines.end()) 
				     // if no constraint is to be handled
    for (unsigned int row=0; row!=n_rows; ++row)
      new_line.push_back (row);
  else
    for (unsigned int row=0; row!=n_rows; ++row)
      if (row == next_constraint->line)
	{
					   // this line is constrained
	  new_line.push_back (-1);
					   // note that @p lines is ordered
	  ++shift;
	  ++next_constraint;
	  if (next_constraint == lines.end())
					     // nothing more to do; finish rest
					     // of loop
	    {
	      for (unsigned int i=row+1; i<n_rows; ++i)
		new_line.push_back (i-shift);
	      break;
	    };
	}
      else
	new_line.push_back (row-shift);


  next_constraint = lines.begin();
				   // note: in this loop we need not check
				   // whether @p next_constraint is a valid
				   // iterator, since @p next_constraint is
				   // only evaluated so often as there are
				   // entries in new_line[*] which tells us
				   // which constraints exist
  for (unsigned int row=0; row<uncondensed.size(); ++row)
    if (new_line[row] != -1)
				       // line not constrained
				       // copy entry
      condensed(new_line[row]) += uncondensed(row);

    else
				       // line must be distributed
      {
	for (unsigned int q=0; q!=next_constraint->entries.size(); ++q) 
	  condensed(new_line[next_constraint->entries[q].first])
	    +=
	    uncondensed(row) * next_constraint->entries[q].second;

	++next_constraint;
      };
}



template <class VectorType>
void
ConstraintMatrix::condense (VectorType &vec) const
{
  Assert (sorted == true, ExcMatrixNotClosed());

                                   // distribute all entries, and set them to zero
  std::vector<ConstraintLine>::const_iterator constraint_line = lines.begin();
  for (; constraint_line!=lines.end(); ++constraint_line)
    {
      for (unsigned int q=0; q!=constraint_line->entries.size(); ++q) 
        vec(constraint_line->entries[q].first)
          += (vec(constraint_line->line) *
              constraint_line->entries[q].second);
      vec(constraint_line->line) = 0.;

				   // in case the constraint is
				   // inhomogeneous, this function is not
				   // appropriate. Throw an exception.
      Assert (constraint_line->inhomogeneity == 0.,
	      ExcMessage ("Inhomogeneous constraint cannot be condensed "
			  "without any matrix specified."));
    }
}



template<typename number, class VectorType>
void
ConstraintMatrix::condense (const SparseMatrix<number> &uncondensed,
			    const VectorType           &uncondensed_vector,
			    SparseMatrix<number>       &condensed,
			    VectorType                 &condensed_vector) const
{
				   // check whether we work on real vectors
				   // or we just used a dummy when calling
				   // the other function above.
  const bool use_vectors = (uncondensed_vector.size() == 0 && 
			    condensed_vector.size() == 0) ? false : true;

  const SparsityPattern &uncondensed_struct = uncondensed.get_sparsity_pattern ();
  
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (uncondensed_struct.is_compressed() == true, ExcMatrixNotClosed());
  Assert (condensed.get_sparsity_pattern().is_compressed() == true, ExcMatrixNotClosed());
  Assert (uncondensed_struct.n_rows() == uncondensed_struct.n_cols(),
	  ExcNotQuadratic());
  Assert (condensed.n() == condensed.m(),
	  ExcNotQuadratic());
  Assert (condensed.n()+n_constraints() == uncondensed.n(),
	  ExcDimensionMismatch(condensed.n()+n_constraints(), uncondensed.n()));
  if (use_vectors == true)
    {
      Assert (condensed_vector.size()+n_constraints() == uncondensed_vector.size(),
	      ExcDimensionMismatch(condensed_vector.size()+n_constraints(),
				   uncondensed_vector.size()));
      Assert (condensed_vector.size() == condensed.m(),
	      ExcDimensionMismatch(condensed_vector.size(), condensed.m()));
    }

				   // store for each line of the matrix
				   // its new line number
				   // after compression. If the shift is
				   // -1, this line will be condensed away
  std::vector<int> new_line;

  new_line.reserve (uncondensed_struct.n_rows());

  std::vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                                shift           = 0;
  const unsigned int n_rows = uncondensed_struct.n_rows();

  if (next_constraint == lines.end()) 
				     // if no constraint is to be handled
    for (unsigned int row=0; row!=n_rows; ++row)
      new_line.push_back (row);
  else
    for (unsigned int row=0; row!=n_rows; ++row)
      if (row == next_constraint->line)
	{
					   // this line is constrained
	  new_line.push_back (-1);
					   // note that @p lines is ordered
	  ++shift;
	  ++next_constraint;
	  if (next_constraint == lines.end())
					     // nothing more to do; finish rest
					     // of loop
	    {
	      for (unsigned int i=row+1; i<n_rows; ++i)
		new_line.push_back (i-shift);
	      break;
	    };
	}
      else
	new_line.push_back (row-shift);


  next_constraint = lines.begin();

				   // note: in this loop we need not check
				   // whether @p next_constraint is a valid
				   // iterator, since @p next_constraint is
				   // only evaluated so often as there are
				   // entries in new_line[*] which tells us
				   // which constraints exist
  for (unsigned int row=0; row<uncondensed_struct.n_rows(); ++row)
    if (new_line[row] != -1)
      {
				       // line not constrained
				       // copy entries if column will not
				       // be condensed away, distribute
				       // otherwise
	for (unsigned int j=uncondensed_struct.get_rowstart_indices()[row];
	     j<uncondensed_struct.get_rowstart_indices()[row+1]; ++j)
	  if (new_line[uncondensed_struct.get_column_numbers()[j]] != -1)
	    condensed.add (new_line[row], new_line[uncondensed_struct.get_column_numbers()[j]],
			   uncondensed.global_entry(j));
	  else 
	    {
					     // let c point to the
					     // constraint of this column
	      std::vector<ConstraintLine>::const_iterator c = lines.begin();
	      while (c->line != uncondensed_struct.get_column_numbers()[j])
		++c;

	      for (unsigned int q=0; q!=c->entries.size(); ++q)
					       // distribute to rows with
					       // appropriate weight
		condensed.add (new_line[row], new_line[c->entries[q].first],
			       uncondensed.global_entry(j) * c->entries[q].second);

				   // take care of inhomogeneity:
				   // need to subtract this element from the
				   // vector. this corresponds to an
				   // explicit elimination in the respective
				   // row of the inhomogeneous constraint in
				   // the matrix with Gauss elimination
	      if (use_vectors == true)
		condensed_vector(new_line[row]) -= uncondensed.global_entry(j) * 
		                                   c->inhomogeneity;
	    }

	if (use_vectors == true)
	  condensed_vector(new_line[row]) += uncondensed_vector(row);	  
      }
    else
				       // line must be distributed
      {
	for (unsigned int j=uncondensed_struct.get_rowstart_indices()[row];
	     j<uncondensed_struct.get_rowstart_indices()[row+1]; ++j)
					   // for each column: distribute
	  if (new_line[uncondensed_struct.get_column_numbers()[j]] != -1)
					     // column is not constrained
	    for (unsigned int q=0; q!=next_constraint->entries.size(); ++q) 
	      condensed.add (new_line[next_constraint->entries[q].first],
			     new_line[uncondensed_struct.get_column_numbers()[j]],
			     uncondensed.global_entry(j) *
			     next_constraint->entries[q].second);
	
	  else
					     // not only this line but
					     // also this col is constrained
	    {
					       // let c point to the constraint
					       // of this column
	      std::vector<ConstraintLine>::const_iterator c = lines.begin();
	      while (c->line != uncondensed_struct.get_column_numbers()[j])
		++c;
	      
	      for (unsigned int p=0; p!=c->entries.size(); ++p)
		for (unsigned int q=0; q!=next_constraint->entries.size(); ++q)
		  condensed.add (new_line[next_constraint->entries[q].first],
				 new_line[c->entries[p].first],
				 uncondensed.global_entry(j) *
				 next_constraint->entries[q].second *
				 c->entries[p].second);

	      if (use_vectors == true)
		for (unsigned int q=0; q!=next_constraint->entries.size(); ++q)
		  condensed_vector (new_line[next_constraint->entries[q].first])
		    -= uncondensed.global_entry(j) * 
		       next_constraint->entries[q].second *
		       c->inhomogeneity;
	    };

				   // condense the vector
	if (use_vectors == true)
	  for (unsigned int q=0; q!=next_constraint->entries.size(); ++q) 
	    condensed_vector(new_line[next_constraint->entries[q].first])
	      +=
	      uncondensed_vector(row) * next_constraint->entries[q].second;

	++next_constraint;
      };
}



template<typename number, class VectorType>
void
ConstraintMatrix::condense (SparseMatrix<number> &uncondensed,
			    VectorType           &vec) const
{
				   // check whether we work on real vectors
				   // or we just used a dummy when calling
				   // the other function above.
  const bool use_vectors = vec.size() == 0 ? false : true;

  const SparsityPattern &sparsity = uncondensed.get_sparsity_pattern ();

  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());
  if (use_vectors == true)
    {
      Assert (vec.size() == sparsity.n_rows(), 
	      ExcDimensionMismatch(vec.size(), sparsity.n_rows()));
    }

  double average_diagonal = 0;
  for (unsigned int i=0; i<uncondensed.m(); ++i)
    average_diagonal += std::fabs (uncondensed.diag_element(i));
  average_diagonal /= uncondensed.m();
  
				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // invalid_unsigned_int, no distribution is
				   // necessary.  otherwise, the number states
				   // which line in the constraint matrix
				   // handles this index
  std::vector<unsigned int> distribute (sparsity.n_rows(),
                                        numbers::invalid_unsigned_int);
  
  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_unsigned_int)
					 // regular line. loop over cols
        {
          for (typename SparseMatrix<number>::iterator
                 entry = uncondensed.begin(row);
               entry != uncondensed.end(row); ++entry)
            {
              const unsigned int column = entry->column();
              
                                               // end of row reached?
                                               // this should not
                                               // happen, since we only
                                               // operate on compressed
                                               // matrices!
              Assert (column != SparsityPattern::invalid_entry,
                      ExcMatrixNotClosed());
	    
              if (distribute[column] != numbers::invalid_unsigned_int)
                                                 // distribute entry at
                                                 // regular row @p row
                                                 // and irregular column
                                                 // sparsity.get_column_numbers()[j];
                                                 // set old entry to
                                                 // zero
                {
                  for (unsigned int q=0;
                       q!=lines[distribute[column]].entries.size(); ++q)
                    uncondensed.add (row,
                                     lines[distribute[column]].entries[q].first,
                                     entry->value() *
                                     lines[distribute[column]].entries[q].second);

				   // need to subtract this element from the
				   // vector. this corresponds to an
				   // explicit elimination in the respective
				   // row of the inhomogeneous constraint in
				   // the matrix with Gauss elimination
		  if (use_vectors == true)
		    vec(row) -= 
		      entry->value() * lines[distribute[column]].inhomogeneity;

                                                   // set old value to zero
                  entry->value() = 0.;
                }
            }
        }
      else
					 // row must be distributed
        {
          for (typename SparseMatrix<number>::iterator
                 entry = uncondensed.begin(row);
               entry != uncondensed.end(row); ++entry)
            {
              const unsigned int column = entry->column();

                                               // end of row reached?
                                               // this should not
                                               // happen, since we only
                                               // operate on compressed
                                               // matrices!
              Assert (column != SparsityPattern::invalid_entry,
                      ExcMatrixNotClosed());

              if (distribute[column] == numbers::invalid_unsigned_int)
                                                 // distribute entry at
                                                 // irregular row
                                                 // @p row and regular
                                                 // column
                                                 // column. set
                                                 // old entry to zero
                {
                  for (unsigned int q=0;
                       q!=lines[distribute[row]].entries.size(); ++q) 
                    uncondensed.add (lines[distribute[row]].entries[q].first,
                                     column,
                                     entry->value() *
                                     lines[distribute[row]].entries[q].second);

                                                   // set old entry to zero
                  entry->value() = 0.;
                }
              else
                                                 // distribute entry at
                                                 // irregular row @p row and
                                                 // irregular column
                                                 // @p column set old entry
                                                 // to one on main
                                                 // diagonal, zero otherwise
                {
                  for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
		    {
		      for (unsigned int q=0;
			   q!=lines[distribute[column]].entries.size(); ++q)
			uncondensed.add (lines[distribute[row]].entries[p].first,
					 lines[distribute[column]].entries[q].first,
					 entry->value() *
					 lines[distribute[row]].entries[p].second *
					 lines[distribute[column]].entries[q].second);

		      if (use_vectors == true)
			vec(lines[distribute[row]].entries[p].first) -= 
			  entry->value() * lines[distribute[row]].entries[p].second *
			  lines[distribute[column]].inhomogeneity;
		    }

                                                   // set old entry to correct
                                                   // value
                  entry->value() = (row == column ? average_diagonal : 0. );
                }
            }

				   // take care of vector
	  if (use_vectors == true)
	    {
	      for (unsigned int q=0; q!=lines[distribute[row]].entries.size(); ++q) 
		vec(lines[distribute[row]].entries[q].first)
		  += (vec(row) * lines[distribute[row]].entries[q].second);

	      vec(lines[distribute[row]].line) = 0.;
	    }
        }
    }
}



template <typename number, class BlockVectorType>
void
ConstraintMatrix::condense (BlockSparseMatrix<number> &uncondensed,
			    BlockVectorType           &vec) const
{
				   // check whether we work on real vectors
				   // or we just used a dummy when calling
				   // the other function above.
  const bool use_vectors = vec.n_blocks() == 0 ? false : true;

  const unsigned int blocks = uncondensed.n_block_rows();
  
  const BlockSparsityPattern &
    sparsity = uncondensed.get_sparsity_pattern ();

  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.get_column_indices() == sparsity.get_row_indices(),
	  ExcNotQuadratic());

  if (use_vectors == true)
    {
      Assert (vec.size() == sparsity.n_rows(), 
	      ExcDimensionMismatch(vec.size(), sparsity.n_rows()));
      Assert (vec.n_blocks() == sparsity.n_block_rows(),
	      ExcDimensionMismatch(vec.n_blocks(), sparsity.n_block_rows()));
    }

  double average_diagonal = 0;
  for (unsigned int b=0; b<uncondensed.n_block_rows(); ++b)
    for (unsigned int i=0; i<uncondensed.block(b,b).m(); ++i)
      average_diagonal += std::fabs (uncondensed.block(b,b).diag_element(i));
  average_diagonal /= uncondensed.m();

  const BlockIndices &
    index_mapping = sparsity.get_column_indices();
  
				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // numbers::invalid_unsigned_int,
				   // no distribution is necessary.
				   // otherwise, the number states which line
				   // in the constraint matrix handles this
				   // index
  std::vector<unsigned int> distribute (sparsity.n_rows(),
                                        numbers::invalid_unsigned_int);
  
  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
				       // get index of this row
				       // within the blocks
      const std::pair<unsigned int,unsigned int>
	block_index = index_mapping.global_to_local(row);
      const unsigned int block_row = block_index.first;
      
      if (distribute[row] == numbers::invalid_unsigned_int)
					 // regular line. loop over
					 // all columns and see
					 // whether this column must
					 // be distributed
	{

					   // to loop over all entries
					   // in this row, we have to
					   // loop over all blocks in
					   // this blockrow and the
					   // corresponding row
					   // therein
	  for (unsigned int block_col=0; block_col<blocks; ++block_col)
	    {
              for (typename SparseMatrix<number>::iterator
                     entry = uncondensed.block(block_row, block_col).begin(block_index.second);
                   entry != uncondensed.block(block_row, block_col).end(block_index.second);
                   ++entry)
                {
                  const unsigned int global_col
                    = index_mapping.local_to_global(block_col,entry->column());
		    
                  if (distribute[global_col] != numbers::invalid_unsigned_int)
                                                     // distribute entry at
                                                     // regular row @p row
                                                     // and irregular column
                                                     // global_col; set old
                                                     // entry to zero
                    {
                      const double old_value = entry->value ();
			
                      for (unsigned int q=0;
                           q!=lines[distribute[global_col]].entries.size(); ++q)
                        uncondensed.add (row,
                                         lines[distribute[global_col]].entries[q].first,
                                         old_value *
                                         lines[distribute[global_col]].entries[q].second);

				   // need to subtract this element from the
				   // vector. this corresponds to an
				   // explicit elimination in the respective
				   // row of the inhomogeneous constraint in
				   // the matrix with Gauss elimination
		      if (use_vectors == true)
			vec(row) -= entry->value() * 
			            lines[distribute[global_col]].inhomogeneity;

                      entry->value() = 0.;
                    }
                }
	    }
	}
      else
	{
					   // row must be
					   // distributed. split the
					   // whole row into the
					   // chunks defined by the
					   // blocks
	  for (unsigned int block_col=0; block_col<blocks; ++block_col)
	    {
              for (typename SparseMatrix<number>::iterator
                     entry = uncondensed.block(block_row, block_col).begin(block_index.second);
                   entry != uncondensed.block(block_row, block_col).end(block_index.second);
                   ++entry)
                {
                  const unsigned int global_col
                    = index_mapping.local_to_global (block_col, entry->column());
		    
                  if (distribute[global_col] ==
                      numbers::invalid_unsigned_int)
                                                     // distribute
                                                     // entry at
                                                     // irregular
                                                     // row @p row
                                                     // and regular
                                                     // column
                                                     // global_col. set
                                                     // old entry to
                                                     // zero
                    {
                      const double old_value = entry->value();
			  
                      for (unsigned int q=0;
                           q!=lines[distribute[row]].entries.size(); ++q) 
                        uncondensed.add (lines[distribute[row]].entries[q].first,
                                         global_col,
                                         old_value *
                                         lines[distribute[row]].entries[q].second);

                      entry->value() = 0.;
                    }
                  else
                                                     // distribute entry at
                                                     // irregular row @p row
                                                     // and irregular column
                                                     // @p global_col set old
                                                     // entry to one if on
                                                     // main diagonal, zero
                                                     // otherwise
                    {
                      const double old_value = entry->value ();
			  
                      for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
			{
			  for (unsigned int q=0; q!=lines[distribute[global_col]].entries.size(); ++q)
			    uncondensed.add (lines[distribute[row]].entries[p].first,
					     lines[distribute[global_col]].entries[q].first,
					     old_value *
					     lines[distribute[row]].entries[p].second *
					     lines[distribute[global_col]].entries[q].second);

			  if (use_vectors == true)
			    vec(lines[distribute[row]].entries[p].first) -= 
			      old_value * lines[distribute[row]].entries[p].second *
			      lines[distribute[global_col]].inhomogeneity;
			}

                      entry->value() = (row == global_col ? average_diagonal : 0. );
                    }
                }
	    }

	  				   // take care of vector
	  if (use_vectors == true)
	    {
	      for (unsigned int q=0; q!=lines[distribute[row]].entries.size(); ++q) 
		vec(lines[distribute[row]].entries[q].first)
		  += (vec(row) * lines[distribute[row]].entries[q].second);

	      vec(lines[distribute[row]].line) = 0.;
	    }
	}
    }
}



template <class VectorType>
void
ConstraintMatrix::set_zero (VectorType &vec) const
{
  Assert (sorted == true, ExcMatrixNotClosed());

  std::vector<ConstraintLine>::const_iterator constraint_line = lines.begin();
  for (; constraint_line!=lines.end(); ++constraint_line)
    vec(constraint_line->line) = 0.;
}



template <typename VectorType>
void
ConstraintMatrix::
distribute_local_to_global (const Vector<double>            &local_vector,
                            const std::vector<unsigned int> &local_dof_indices,
                            VectorType                      &global_vector) const
{
  Assert (local_vector.size() == local_dof_indices.size(),
          ExcDimensionMismatch(local_vector.size(), local_dof_indices.size()));
  Assert (sorted == true, ExcMatrixNotClosed());

  const unsigned int n_local_dofs = local_vector.size();
  
                                   // have a special case where there are no
                                   // constraints at all, since then we can be
                                   // a lot faster
  Threads::ThreadMutex::ScopedLock lock(mutex);
  for (unsigned int i=0; i<n_local_dofs; ++i)
    {
				   // let's see if we can use the bool
				   // vector that tells about whether a
				   // certain constraint exists.
      if (constraint_line_exists.size() <= local_dof_indices[i])
	{
	  global_vector(local_dof_indices[i]) += local_vector(i);
	  continue;
	}
      if (constraint_line_exists[local_dof_indices[i]] == false)
	{
	  global_vector(local_dof_indices[i]) += local_vector(i);
	  continue;
	}
                                           // first figure out whether this
                                           // dof is constrained
      ConstraintLine index_comparison;
      index_comparison.line = local_dof_indices[i];

      const std::vector<ConstraintLine>::const_iterator
	position = std::lower_bound (lines.begin(),
				     lines.end(),
				     index_comparison);

                                           // if the line is not
                                           // constrained, then simply
                                           // copy the data. otherwise
                                           // distribute it, but make
                                           // sure we don't touch the
                                           // entries of fixed dofs
					   //
					   // there is one critical
					   // point: sometimes a dof
					   // may be both constrained
					   // and fixed, for example
					   // hanging nodes in 3d at
					   // the boundary. in that
					   // case, we don't quite
					   // know what to do --
					   // handle the constraint or
					   // the fixed
					   // value. however, this
					   // isn't so hard if all the
					   // nodes that this node is
					   // constrained to are also
					   // fixed nodes, in which
					   // case we could do both
					   // but opt to copy the
					   // element. however, we
					   // have to check that all
					   // the nodes to which it is
					   // constrained are also
					   // fixed
      Assert (position->line == local_dof_indices[i],
	      ExcInternalError());
      for (unsigned int j=0; j<position->entries.size(); ++j)
	global_vector(position->entries[j].first)
	  += local_vector(i) * position->entries[j].second;
    }
}



template <typename MatrixType>
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<double>        &local_matrix,
                            const std::vector<unsigned int> &local_dof_indices,
                            MatrixType                      &global_matrix) const
{
				   // create a dummy and hand on to the
				   // function actually implementing this
				   // feature further down.
  Vector<double> dummy(0);
  distribute_local_to_global (local_matrix, dummy, local_dof_indices,
			      global_matrix, dummy, 
			      internal::bool2type<IsBlockMatrix<MatrixType>::value>());
}



template <typename MatrixType, typename VectorType>
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<double>        &local_matrix,
			    const Vector<double>            &local_vector,
                            const std::vector<unsigned int> &local_dof_indices,
                            MatrixType                      &global_matrix,
			    VectorType                      &global_vector) const
{
				   // enter the internal function with the
				   // respective block information set, the
				   // actual implementation follows further
				   // down.
  distribute_local_to_global (local_matrix, local_vector, local_dof_indices,
			      global_matrix, global_vector, 
			      internal::bool2type<IsBlockMatrix<MatrixType>::value>());
}



template <typename SparsityType>
void
ConstraintMatrix::
add_entries_local_to_global (const std::vector<unsigned int> &local_dof_indices,
			     SparsityType                    &sparsity_pattern,
			     const bool                       keep_constrained_entries,
			     const Table<2,bool>             &dof_mask) const
{
				   // enter the internal function with the
				   // respective block information set, the
				   // actual implementation follows further
				   // down.
  add_entries_local_to_global (local_dof_indices, sparsity_pattern,
			       keep_constrained_entries, dof_mask,
			       internal::bool2type<IsBlockMatrix<SparsityType>::value>());
}



template<class VectorType>
void
ConstraintMatrix::distribute (const VectorType &condensed,
			      VectorType       &uncondensed) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (condensed.size()+n_constraints() == uncondensed.size(),
	  ExcDimensionMismatch(condensed.size()+n_constraints(),
			       uncondensed.size()));

				   // store for each line of the new vector
				   // its old line number before
				   // distribution. If the shift is
				   // -1, this line was condensed away
  std::vector<int> old_line;

  old_line.reserve (uncondensed.size());

  std::vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                                shift           = 0;
  unsigned int n_rows = uncondensed.size();

  if (next_constraint == lines.end()) 
				     // if no constraint is to be handled
    for (unsigned int row=0; row!=n_rows; ++row)
      old_line.push_back (row);
  else
    for (unsigned int row=0; row!=n_rows; ++row)
      if (row == next_constraint->line)
	{
					   // this line is constrained
	  old_line.push_back (-1);
					   // note that @p lines is ordered
	  ++shift;
	  ++next_constraint;
	  if (next_constraint == lines.end())
					     // nothing more to do; finish rest
					     // of loop
	    {
	      for (unsigned int i=row+1; i<n_rows; ++i)
		old_line.push_back (i-shift);
	      break;
	    };
	}
      else
	old_line.push_back (row-shift);


  next_constraint = lines.begin();
				   // note: in this loop we need not check
				   // whether @p next_constraint is a valid
				   // iterator, since @p next_constraint is
				   // only evaluated so often as there are
				   // entries in new_line[*] which tells us
				   // which constraints exist
  for (unsigned int line=0; line<uncondensed.size(); ++line) 
    if (old_line[line] != -1)
				       // line was not condensed away
      uncondensed(line) = condensed(old_line[line]);
    else
      {
					 // line was condensed away,
					 // create it newly. first set
					 // it to zero
	uncondensed(line) = next_constraint->inhomogeneity;
					 // then add the different
					 // contributions
	for (unsigned int i=0; i<next_constraint->entries.size(); ++i)
	  uncondensed(line) += (condensed(old_line[next_constraint->entries[i].first]) *
				next_constraint->entries[i].second);
	++next_constraint;
      };
}



template<class VectorType>
void
ConstraintMatrix::distribute (VectorType &vec) const
{
  Assert (sorted == true, ExcMatrixNotClosed());

  std::vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  for (; next_constraint != lines.end(); ++next_constraint) 
    {
				       // fill entry in line
				       // next_constraint.line by adding the
				       // different contributions
      double new_value = next_constraint->inhomogeneity;
      for (unsigned int i=0; i<next_constraint->entries.size(); ++i)
	new_value += (vec(next_constraint->entries[i].first) *
                      next_constraint->entries[i].second);
      vec(next_constraint->line) = new_value;
    }
}


				   // Some helper definitions for the
				   // local_to_global functions.
namespace internals
{
				   // this struct contains all the
				   // information we need to store about
				   // which global entries (global_row) are
				   // given rise by local entries
				   // (local_row) or some constraints.
  struct distributing 
  {
    distributing (const unsigned int global_row = deal_II_numbers::invalid_unsigned_int,
		  const unsigned int local_row = deal_II_numbers::invalid_unsigned_int);
    distributing (const distributing &in);
    ~distributing ();
    distributing & operator = (const distributing &in);
    bool operator < (const distributing &in);
    unsigned int global_row;
    unsigned int local_row;
    mutable std::vector<std::pair<unsigned int,double> > *constraints;
  };

  inline
  distributing::distributing (const unsigned int global_row,
			      const unsigned int local_row) :
    global_row (global_row),
    local_row (local_row),
    constraints (0) {}

  inline
  distributing::distributing (const distributing &in) :
    constraints (0)
  {*this = (in);}

  inline
  distributing::~distributing ()
  {
    if (constraints != 0)
      {
	delete constraints;
	constraints = 0;
      }
  }

  inline
  distributing & distributing::operator = (const distributing &in)
  {
    global_row = in.global_row;
    local_row = in.local_row;
				   // the constraints pointer should not
				   // contain any data here.
    Assert (constraints == 0, ExcInternalError());

    if (in.constraints != 0)
      {
	constraints = in.constraints;
	in.constraints = 0;
      }
    return *this;
  }

  inline
  bool distributing::operator < (const distributing &in)
  {
    return global_row < in.global_row;
  }

				   // a function that appends an additional
				   // row to the list of values, or appends
				   // a value to an already existing
				   // row. Similar functionality as for
				   // std::map<unsigned int,distributing>,
				   // but here done for a std::vector of
				   // data type distributing, and much
				   // faster.
  inline
  void
  insert_index (std::vector<distributing> &my_indices,
		const unsigned int row,
		const std::pair<unsigned int,double> constraint)
  {
    typedef std::vector<distributing>::iterator index_iterator;
    index_iterator pos, pos1;
    distributing row_value (row);

				   // check whether the list was really
				   // sorted before entering here
#ifdef DEBUG
    for (unsigned int i=1; i<my_indices.size(); ++i)
      Assert (my_indices[i-1] < my_indices[i], ExcInternalError());
#endif

    if (my_indices.size() == 0 || my_indices.back().global_row < row)
      {
	my_indices.push_back(row_value);
	pos1 = my_indices.end()-1;
      }
    else
      {
	pos = std::lower_bound (my_indices.begin(),my_indices.end(), row_value);
	if (pos->global_row == row)
	  pos1 = pos;
	else
	  pos1 = my_indices.insert(pos, row_value);
      }

    if (&*pos1->constraints == 0)
      pos1->constraints = 
	new std::vector<std::pair<unsigned int,double> > (1,constraint);
    else
      pos1->constraints->push_back (constraint);
  }


				   // this sort algorithm sorts a vector of
				   // distributing elements, but does not
				   // take the constraints into
				   // account. this means that in case that
				   // constraints are already inserted, this
				   // function does not work as
				   // expected. shellsort is very fast in
				   // case the indices are already sorted
				   // (which is the usual case with DG
				   // elements), and not too slow in other
				   // cases
  inline
  void
  list_shellsort (std::vector<distributing> &my_indices)
  {
    unsigned int i, j, j2, temp, templ, istep;
    unsigned step;

				   // in debug mode, check whether the
				   // constraints are really empty.
#ifdef DEBUG
    for (unsigned int i=0; i<my_indices.size(); ++i)
      Assert (&*my_indices[i].constraints == 0, ExcInternalError());
#endif

    const unsigned int length = my_indices.size();
    step = length/2;
    while (step > 0)
      {
	for (i=step; i < length; i++)
	  {
	    istep = step;
	    j = i;
	    j2 = j-istep;
	    temp = my_indices[i].global_row;
	    templ = my_indices[i].local_row;
	    if (my_indices[j2].global_row > temp) 
	      {
		while ((j >= istep) && (my_indices[j2].global_row > temp))
		  {
		    my_indices[j].global_row = my_indices[j2].global_row;
		    my_indices[j].local_row = my_indices[j2].local_row;
		    j = j2;
		    j2 -= istep;
		  }
		my_indices[j].global_row = temp;
		my_indices[j].local_row = templ;
	      }
	  }
	step = step>>1;
      }
  }

				   // this function for block matrices
				   // creates some block indices for the
				   // list of local dofs and transforms the
				   // indices to local block indices.
  template <class BlockType>
  inline
  void
  make_block_starts (const BlockType           &block_object,
		     std::vector<unsigned int> &row_indices,
		     std::vector<unsigned int> &block_starts)
  {
    Assert (block_starts.size() == block_object.n_block_rows() + 1,
	    ExcDimensionMismatch(block_starts.size(),
				 block_object.n_block_rows()+1));

    typedef std::vector<unsigned int>::iterator row_iterator;
    row_iterator col_indices = row_indices.begin();

    const unsigned int num_blocks = block_object.n_block_rows();

				   // find end of rows.
    block_starts[0] = 0;
    for (unsigned int i=1;i<num_blocks;++i)
      {
	row_iterator first_block = 
	  std::lower_bound (col_indices,
			    row_indices.end(),
			    block_object.get_row_indices().block_start(i));
	block_starts[i] = first_block - row_indices.begin();
	col_indices = first_block;
      }

				   // transform row indices to local index
				   // space
    for (unsigned int i=block_starts[1]; i<row_indices.size(); ++i)
      row_indices[i] = block_object.get_row_indices().global_to_local(row_indices[i]).second;
  }

}


				   // internal implementation for
				   // distribute_local_to_global for
				   // standard (non-block) matrices
template <typename MatrixType, typename VectorType>
inline
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<double>        &local_matrix,
			    const Vector<double>            &local_vector,
                            const std::vector<unsigned int> &local_dof_indices,
                            MatrixType                      &global_matrix,
			    VectorType                      &global_vector,
			    internal::bool2type<false>) const
{
				   // check whether we work on real vectors
				   // or we just used a dummy when calling
				   // the other function above.
  const bool use_vectors = (local_vector.size() == 0 && 
			    global_vector.size() == 0) ? false : true;

  Assert (local_matrix.n() == local_dof_indices.size(),
          ExcDimensionMismatch(local_matrix.n(), local_dof_indices.size()));
  Assert (local_matrix.m() == local_dof_indices.size(),
          ExcDimensionMismatch(local_matrix.m(), local_dof_indices.size()));
  Assert (global_matrix.m() == global_matrix.n(), ExcNotQuadratic());
  if (use_vectors == true)
    {
      Assert (local_matrix.m() == local_vector.size(),
	      ExcDimensionMismatch(local_matrix.m(), local_vector.size()));
      Assert (global_matrix.m() == global_vector.size(),
	      ExcDimensionMismatch(global_matrix.m(), global_vector.size()));
    }
  Assert (sorted == true, ExcMatrixNotClosed());

  const unsigned int n_local_dofs = local_dof_indices.size();

  double average_diagonal = 0;
  for (unsigned int i=0; i<n_local_dofs; ++i)
    average_diagonal += std::fabs (local_matrix(i,i));
  average_diagonal /= n_local_dofs;

				   // when distributing the local data to
				   // the global matrix, we can quite
				   // cheaply sort the indices (obviously,
				   // this introduces the need for
				   // allocating some memory on the way, but
				   // we need to do this only for rows,
				   // whereas the distribution process
				   // itself goes over rows and
				   // columns). This has the advantage that
				   // when writing into the global matrix,
				   // we can make use of the sortedness.

				   // so the first step is to create a
				   // sorted list of all row values that are
				   // possible. these values are either the
				   // rows from unconstrained dofs, or some
				   // indices introduced by dofs constrained
				   // to a combination of some other
				   // dofs. regarding the data type, choose
				   // an STL vector of a pair of unsigned
				   // ints (for global columns) and internal
				   // data (containing local columns +
				   // possible jumps from
				   // constraints). Choosing an STL map or
				   // anything else M.K. knows of would be
				   // much more expensive here!
  std::vector<internals::distributing> my_indices (n_local_dofs);
  std::vector<std::pair<unsigned int, const ConstraintLine *> > constraint_lines;

				   // cache whether we have to resolve any
				   // indirect rows generated from resolving
				   // constrained dofs.
  bool have_indirect_rows = false;
  {
    unsigned int added_rows = 0;
				   // first add the indices in an unsorted
				   // way and only keep track of the
				   // constraints that appear. They are
				   // resolved in a second step.
    for (unsigned int i = 0; i<n_local_dofs; ++i)
      {
	if (constraint_line_exists.size() <= local_dof_indices[i] ||
	    constraint_line_exists[local_dof_indices[i]] == false)
	  {
	    my_indices[added_rows].global_row = local_dof_indices[i];
	    my_indices[added_rows].local_row = i;
	    ++added_rows;
	    continue;
	  }

	ConstraintLine index_comparison;
	index_comparison.line = local_dof_indices[i];

	const std::vector<ConstraintLine>::const_iterator
	  position = std::lower_bound (lines.begin(),
				       lines.end(),
				       index_comparison);
	Assert (position->line == local_dof_indices[i],
		ExcInternalError());

	constraint_lines.push_back (std::make_pair<unsigned int,
				    const ConstraintLine *>(i,&*position));
      }
    Assert (constraint_lines.size() + added_rows == n_local_dofs,
	    ExcInternalError());
    my_indices.resize (added_rows);
  }
  internals::list_shellsort (my_indices);

  Threads::ThreadMutex::ScopedLock lock(mutex);

				   // now in the second step actually
				   // resolve the constraints
  const unsigned int n_constrained_dofs = constraint_lines.size();
  for (unsigned int i=0; i<n_constrained_dofs; ++i)
    {
      const unsigned int local_row = constraint_lines[i].first;
      const unsigned int global_row = local_dof_indices[local_row];
      const ConstraintLine * position = constraint_lines[i].second;
      for (unsigned int q=0; q<position->entries.size(); ++q)
	{
	  have_indirect_rows = true;
	  internals::insert_index(my_indices, position->entries[q].first,
				  std::make_pair<unsigned int,double> 
				  (local_row, position->entries[q].second));
	}

				   // to make sure that the global matrix
				   // remains invertible, we need to do
				   // something with the diagonal
				   // elements. add the absolute value of
				   // the local matrix, so the resulting
				   // entry will always be positive and
				   // furthermore be in the same order of
				   // magnitude as the other elements of the
				   // matrix
				   //
				   // note that this also captures the
				   // special case that a dof is both
				   // constrained and fixed (this can happen
				   // for hanging nodes in 3d that also
				   // happen to be on the boundary). in that
				   // case, following the above program
				   // flow, it is realized that when
				   // distributing the row and column no
				   // elements of the matrix are actually
				   // touched if all the degrees of freedom
				   // to which this dof is constrained are
				   // also constrained (the usual case with
				   // hanging nodes in 3d). however, in the
				   // line below, we do actually do
				   // something with this dof
      const double new_diagonal = std::fabs(local_matrix(local_row,local_row)) != 0 ?
	std::fabs(local_matrix(local_row,local_row)) : average_diagonal;
      global_matrix.add(global_row, global_row, new_diagonal);
    }

  const unsigned int n_actual_dofs = my_indices.size();

				   // create arrays for the column data
				   // (indices and values) that will then be
				   // written into the matrix.
  std::vector<unsigned int> cols (n_actual_dofs);
  std::vector<double>       vals (n_actual_dofs);

  typedef std::vector<std::pair<unsigned int,double> > constraint_format;

				   // now do the actual job. 
  for (unsigned int i=0; i<n_actual_dofs; ++i)
    {
      const unsigned int row = my_indices[i].global_row;
      const unsigned int loc_row = my_indices[i].local_row;
      unsigned int * col_ptr = &cols[0];
      double * val_ptr = &vals[0];
      double val = 0;

				   // fast function if there are no indirect
				   // references to any of the local rows at
				   // all on this set of dofs (saves a lot
				   // of checks). the only check we actually
				   // need to perform is whether the matrix
				   // element is zero. 
      if (have_indirect_rows == false)
	{
	  Assert(loc_row < n_local_dofs, ExcInternalError());
	  const double * matrix_ptr = &local_matrix(loc_row, 0);

	  for (unsigned int j=0; j < n_actual_dofs; ++j)
	    {
	      const unsigned int loc_col = my_indices[j].local_row;
	      Assert(loc_col < n_local_dofs, ExcInternalError());

	      const double col_val = matrix_ptr[loc_col];
	      if (col_val != 0)
		{
		  *val_ptr++ = col_val;
		  *col_ptr++ = my_indices[j].global_row;
		}
	    }

	  if (use_vectors == true)
	    {
	      val = local_vector(loc_row);

				   // need to account for inhomogeneities
				   // here: thie corresponds to eliminating
				   // the respective column in the local
				   // matrix with value on the right hand
				   // side.
	      for (unsigned int i=0; i<constraint_lines.size(); ++i)
		val -= constraint_lines[i].second->inhomogeneity *
		       matrix_ptr[constraint_lines[i].first];
	    }
	}

				   // more difficult part when there are
				   // indirect references and when we need
				   // to do some more checks.
      else
	{
	  const double * matrix_ptr = 0;
	  if (loc_row != deal_II_numbers::invalid_unsigned_int)
	    {
	      Assert (loc_row < n_local_dofs, ExcInternalError());
	      matrix_ptr = &local_matrix(loc_row, 0);
	    }
	  for (unsigned int j=0; j < n_actual_dofs; ++j)
	    {
	      double col_val;
	      const unsigned int loc_col = my_indices[j].local_row;

				   // case 1: row has direct contribution in
				   // local matrix
	      if (loc_row != deal_II_numbers::invalid_unsigned_int)
		{
				   // case 1a: col has direct contribution
				   // in local matrix
		  if (loc_col != deal_II_numbers::invalid_unsigned_int)
		    {
		      Assert (loc_col < n_local_dofs, ExcInternalError());
		      col_val = matrix_ptr[loc_col];
		    }
				   // case 1b: col has no direct
				   // contribution in local matrix
		  else
		    col_val = 0;

				   // account for indirect contributions by
				   // constraints
		  if (my_indices[j].constraints != 0)
		    {
		      constraint_format &constraint_j = *my_indices[j].constraints;
		      for (unsigned int p=0; p<constraint_j.size(); ++p)
			col_val += matrix_ptr[constraint_j[p].first] 
		                   *
		                   constraint_j[p].second;
		    }
		}

				   // case 2: row has no direct contribution in
				   // local matrix
	      else
		col_val = 0;

				   // account for indirect contributions by
				   // constraints in row, going trough the
				   // direct and indirect references in the
				   // given column.
	      if (my_indices[i].constraints != 0)
		{
		  constraint_format &constraint_i = *my_indices[i].constraints;
		  Assert (constraint_i.size() > 0, ExcInternalError());

		  for (unsigned int q=0; q<constraint_i.size(); ++q)
		    {
		      double add_this = loc_col != deal_II_numbers::invalid_unsigned_int ?
			local_matrix(constraint_i[q].first, loc_col) : 0;
  
		      if (my_indices[j].constraints != 0)
			{
			  constraint_format &constraint_j = *my_indices[j].constraints;
			  for (unsigned int p=0; p<constraint_j.size(); ++p)
			    add_this += local_matrix(constraint_i[q].first,
						     constraint_j[p].first) 
			                *
			                constraint_j[p].second;
			}
		      col_val += add_this * constraint_i[q].second;
		    }
		}

		   

				   // if we got some nontrivial value,
				   // append it to the array of values.
	      if (col_val != 0)
		{
		  *val_ptr++ = col_val;
		  *col_ptr++ = my_indices[j].global_row;
		}
	    }

				   // now to the vectors. besides doing the
				   // same job as we did above (i.e.,
				   // distribute the content of the local
				   // vector into the global one), need to
				   // account for inhomogeneities here: thie
				   // corresponds to eliminating the
				   // respective column in the local matrix
				   // with value on the right hand side.
	  if (use_vectors == true)
	    {
	      if (loc_row != deal_II_numbers::invalid_unsigned_int)
		{
		  Assert (loc_row < n_local_dofs,
			  ExcInternalError());
		  val = local_vector(loc_row);
		  for (unsigned int i=0; i<constraint_lines.size(); ++i)
		    val -= constraint_lines[i].second->inhomogeneity *
		           matrix_ptr[constraint_lines[i].first];
		}

	      if (my_indices[i].constraints != 0)
		{
		  std::vector<std::pair<unsigned int,double> > &constraint_i = 
		    *my_indices[i].constraints;

		  for (unsigned int q=0; q<constraint_i.size(); ++q)
		    {
		      const unsigned int loc_row_q = constraint_i[q].first;
		      double add_this = local_vector (loc_row_q);
		      for (unsigned int k=0; k<constraint_lines.size(); ++k)
			add_this -= constraint_lines[k].second->inhomogeneity *
		                    local_matrix(loc_row_q,constraint_lines[k].first);
		      val += add_this * constraint_i[q].second;
		    }
		}
	    }
	}

				   // finally, write all the information
				   // that accumulated under the given
				   // process into the global matrix row and
				   // into the vector
      const unsigned int n_values = col_ptr - &cols[0];
      Assert (n_values == (unsigned int)(val_ptr - &vals[0]),
	      ExcInternalError());
      if (n_values > 0)
	global_matrix.add(row, n_values, &cols[0], &vals[0], false, true);
      if (val != 0)
	global_vector(row) += val;
    }
}



template <typename MatrixType, typename VectorType>
inline
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<double>        &local_matrix,
			    const Vector<double>            &local_vector,
                            const std::vector<unsigned int> &local_dof_indices,
                            MatrixType                      &global_matrix,
			    VectorType                      &global_vector,
			    internal::bool2type<true>) const
{
				   // similar function as above, but now
				   // specialized for block matrices. See
				   // the other function for additional
				   // comments.

  const bool use_vectors = (local_vector.size() == 0 && 
			    global_vector.size() == 0) ? false : true;

  Assert (local_matrix.n() == local_dof_indices.size(),
          ExcDimensionMismatch(local_matrix.n(), local_dof_indices.size()));
  Assert (local_matrix.m() == local_dof_indices.size(),
          ExcDimensionMismatch(local_matrix.m(), local_dof_indices.size()));
  Assert (global_matrix.m() == global_matrix.n(), ExcNotQuadratic());
  Assert (global_matrix.n_block_rows() == global_matrix.n_block_cols(),
	  ExcNotQuadratic());
  if (use_vectors == true)
    {
      Assert (local_matrix.m() == local_vector.size(),
	      ExcDimensionMismatch(local_matrix.m(), local_vector.size()));
      Assert (global_matrix.m() == global_vector.size(),
	      ExcDimensionMismatch(global_matrix.m(), global_vector.size()));
    }
  Assert (sorted == true, ExcMatrixNotClosed());

  const unsigned int n_local_dofs = local_dof_indices.size();
  const unsigned int num_blocks   = global_matrix.n_block_rows();

  double average_diagonal = 0;
  for (unsigned int i=0; i<n_local_dofs; ++i)
    average_diagonal += std::fabs (local_matrix(i,i));
  average_diagonal /= n_local_dofs;

  std::vector<internals::distributing> my_indices (n_local_dofs);
  std::vector<std::pair<unsigned int, const ConstraintLine *> > constraint_lines;

  bool have_indirect_rows = false;
  {
    unsigned int added_rows = 0;
    for (unsigned int i = 0; i<n_local_dofs; ++i)
      {
	if (constraint_line_exists.size() <= local_dof_indices[i] ||
	    constraint_line_exists[local_dof_indices[i]] == false)
	  {
	    my_indices[added_rows].global_row = local_dof_indices[i];
	    my_indices[added_rows].local_row = i;
	    ++added_rows;
	    continue;
	  }

	ConstraintLine index_comparison;
	index_comparison.line = local_dof_indices[i];

	const std::vector<ConstraintLine>::const_iterator
	  position = std::lower_bound (lines.begin(),
				       lines.end(),
				       index_comparison);
	Assert (position->line == local_dof_indices[i],
		ExcInternalError());

	constraint_lines.push_back (std::make_pair<unsigned int,
				    const ConstraintLine *>(i,&*position));
      }
    Assert (constraint_lines.size() + added_rows == n_local_dofs,
	    ExcInternalError());
    my_indices.resize (added_rows);
  }
  internals::list_shellsort (my_indices);

  Threads::ThreadMutex::ScopedLock lock(mutex);

  const unsigned int n_constrained_dofs = constraint_lines.size();
  for (unsigned int i=0; i<n_constrained_dofs; ++i)
    {
      const unsigned int local_row = constraint_lines[i].first;
      const unsigned int global_row = local_dof_indices[local_row];
      const ConstraintLine * position = constraint_lines[i].second;
      for (unsigned int q=0; q<position->entries.size(); ++q)
	{
	  have_indirect_rows = true;
	  internals::insert_index(my_indices, position->entries[q].first,
				  std::make_pair<unsigned int,double> 
				  (local_row, position->entries[q].second));
	}

      const double new_diagonal = std::fabs(local_matrix(local_row,local_row)) != 0 ?
	std::fabs(local_matrix(local_row,local_row)) : average_diagonal;
      global_matrix.add(global_row, global_row, new_diagonal);
    }

  const unsigned int n_actual_dofs = my_indices.size();

  std::vector<unsigned int> localized_indices (n_actual_dofs);
  for (unsigned int i=0; i<n_actual_dofs; ++i)
    localized_indices[i] = my_indices[i].global_row;

				   // additional construct that also takes
				   // care of block indices.
  std::vector<unsigned int> block_starts(num_blocks+1, n_actual_dofs);
  internals::make_block_starts (global_matrix, localized_indices, block_starts);

  std::vector<unsigned int> cols (n_actual_dofs);
  std::vector<double>       vals (n_actual_dofs);
  typedef std::vector<std::pair<unsigned int,double> > constraint_format;

				   // the basic difference to the
				   // non-block variant from now onwards
				   // is that we go through the blocks
				   // of the matrix separately.
  for (unsigned int block=0; block<num_blocks; ++block)
    {
      const unsigned int next_block = block_starts[block+1];
      for (unsigned int i=block_starts[block]; i<next_block; ++i)
	{
	  const unsigned int row = localized_indices[i];
	  const unsigned int loc_row = my_indices[i].local_row;

	  for (unsigned int block_col=0; block_col<num_blocks; ++block_col)
	    {
	      const unsigned int next_block_col = block_starts[block_col+1];
	      unsigned int * col_ptr = &cols[0];
	      double * val_ptr = &vals[0];
	      if (have_indirect_rows == false)
		{
		  Assert(loc_row < n_local_dofs, ExcInternalError());
		  const double * matrix_ptr = &local_matrix(loc_row, 0);

		  for (unsigned int j=block_starts[block_col]; j < next_block_col; ++j)
		    {
		      const unsigned int loc_col = my_indices[j].local_row;
		      Assert(loc_col < n_local_dofs, ExcInternalError());

		      const double col_val = matrix_ptr[loc_col];
		      if (col_val != 0)
			{
			  *val_ptr++ = col_val;
			  *col_ptr++ = localized_indices[j];
			}
		    }
		}

	      else
		{
		  const double * matrix_ptr = 0;
		  if (loc_row != deal_II_numbers::invalid_unsigned_int)
		    {
		      Assert (loc_row < n_local_dofs, ExcInternalError());
		      matrix_ptr = &local_matrix(loc_row, 0);
		    }
		  for (unsigned int j=block_starts[block_col]; j < next_block_col; ++j)
		    {
		      double col_val;
		      const unsigned int loc_col = my_indices[j].local_row;

		      if (loc_row != deal_II_numbers::invalid_unsigned_int)
			{
			  col_val = loc_col != deal_II_numbers::invalid_unsigned_int ?
			    matrix_ptr[loc_col] : 0;

				   // account for indirect contributions by
				   // constraints
			  if (my_indices[j].constraints != 0)
			    {
			      constraint_format &constraint_j = 
				*my_indices[j].constraints;

			      for (unsigned int p=0; p<constraint_j.size(); ++p)
				col_val += local_matrix(loc_row,
						        constraint_j[p].first) 
				           *
				           constraint_j[p].second;
			    }
			}

		      else
			col_val = 0;

		      if (my_indices[i].constraints != 0)
			{
			  constraint_format &constraint_i = *my_indices[i].constraints;

			  for (unsigned int q=0; q<constraint_i.size(); ++q)
			    {
			      double add_this = 
				loc_col != deal_II_numbers::invalid_unsigned_int ?
				local_matrix(constraint_i[q].first, loc_col) : 0;
  
			      if (my_indices[j].constraints != 0)
				{
				  constraint_format &constraint_j = 
				    *my_indices[j].constraints;

				  for (unsigned int p=0; p<constraint_j.size(); ++p)
				    add_this += local_matrix(constraint_i[q].first,
							     constraint_j[p].first) 
				                *
				                constraint_j[p].second;
				}
			      col_val += add_this * constraint_i[q].second;
			    }
			}

		      if (col_val != 0)
			{
			  *col_ptr++ = localized_indices[j];
			  *val_ptr++ = col_val;
			}
		    }
		}

				   // finally, write all the information
				   // that accumulated under the given
				   // process into the global matrix row and
				   // into the vector. For the block matrix,
				   // go trough the individual blocks and
				   // look which entries we need to set.
	      const unsigned int n_values = col_ptr - &cols[0];
	      Assert (n_values == (unsigned int)(val_ptr - &vals[0]),
		      ExcInternalError());
	      if (n_values > 0)
		global_matrix.block(block, block_col).add(row, n_values,
							  &cols[0], &vals[0],
							  false, true);
	    }

	  if (use_vectors == true)
	    {
	      double val = 0;
	      if (loc_row != deal_II_numbers::invalid_unsigned_int)
		{
		  Assert (loc_row < n_local_dofs,
			  ExcInternalError());
		  val = local_vector(loc_row);
		  for (unsigned int i=0; i<constraint_lines.size(); ++i)
		    val -= constraint_lines[i].second->inhomogeneity *
		      local_matrix(loc_row,constraint_lines[i].first);
		}

	      if (my_indices[i].constraints != 0)
		{
		  std::vector<std::pair<unsigned int,double> > &constraint_i = 
		    *my_indices[i].constraints;

		  for (unsigned int q=0; q<constraint_i.size(); ++q)
		    {
		      const unsigned int loc_row_q = constraint_i[q].first;
		      double add_this = local_vector (loc_row_q);
		      for (unsigned int k=0; k<constraint_lines.size(); ++k)
			add_this -= constraint_lines[k].second->inhomogeneity *
			            local_matrix(loc_row_q,constraint_lines[k].first);
		      val += add_this * constraint_i[q].second;
		    }
		}
	      if (val != 0)
		global_vector(my_indices[i].global_row) += val;
	    }
	}
    }
}



template <typename SparsityType>
inline
void
ConstraintMatrix::
add_entries_local_to_global (const std::vector<unsigned int> &local_dof_indices,
			     SparsityType                    &sparsity_pattern,
			     const bool                       keep_constrained_entries,
			     const Table<2,bool>             &dof_mask,
			     internal::bool2type<false> ) const
{
  Assert (sparsity_pattern.n_rows() == sparsity_pattern.n_cols(), ExcNotQuadratic());

  const unsigned int n_local_dofs = local_dof_indices.size();
  bool dof_mask_is_active = false;
  if (dof_mask.n_rows() == n_local_dofs)
    {
      dof_mask_is_active = true;
      Assert (dof_mask.n_cols() == n_local_dofs,
	      ExcDimensionMismatch(dof_mask.n_cols(), n_local_dofs));
    }

				   // if the dof mask is not active, all we
				   // have to do is to add some indices in a
				   // matrix format. To do this, we first
				   // create an array of all the indices
				   // that are to be added. these indices
				   // are the local dof indices plus some
				   // indices that come from constraints.
  if (dof_mask_is_active == false)
    {
      std::vector<unsigned int> actual_dof_indices (n_local_dofs);
      unsigned int added_rows = 0;
      bool have_indirect_rows = false;
      std::vector<std::pair<unsigned int, const ConstraintLine *> > constraint_lines;
      for (unsigned int i = 0; i<n_local_dofs; ++i)
	{
	  if (constraint_line_exists.size() <= local_dof_indices[i] ||
	      constraint_line_exists[local_dof_indices[i]] == false)
	    {
	      actual_dof_indices[added_rows] = local_dof_indices[i];
	      ++added_rows;
	      continue;
	    }

	  ConstraintLine index_comparison;
	  index_comparison.line = local_dof_indices[i];

	  const std::vector<ConstraintLine>::const_iterator
	    position = std::lower_bound (lines.begin(),
					 lines.end(),
					 index_comparison);
	  Assert (position->line == local_dof_indices[i],
		  ExcInternalError());

	  constraint_lines.push_back (std::make_pair<unsigned int,
				      const ConstraintLine *>(i,&*position));
      }
      Assert (constraint_lines.size() + added_rows == n_local_dofs,
	      ExcInternalError());
      actual_dof_indices.resize (added_rows);
      std::sort (actual_dof_indices.begin(), actual_dof_indices.end());

      Threads::ThreadMutex::ScopedLock lock(mutex);

      const unsigned int n_constrained_dofs = constraint_lines.size();
      for (unsigned int i=0; i<n_constrained_dofs; ++i)
	{
	  const unsigned int local_row = constraint_lines[i].first;
	  const unsigned int global_row = local_dof_indices[local_row];
	  const ConstraintLine * position = constraint_lines[i].second;
	  for (unsigned int q=0; q<position->entries.size(); ++q)
	    {
	      have_indirect_rows = true;
	      const unsigned int new_index = position->entries[q].first;
	      if (actual_dof_indices.back() < new_index)
		{
		  actual_dof_indices.push_back(new_index);
		}
	      else
		{
		  std::vector<unsigned int>::iterator it = 
		    std::lower_bound(actual_dof_indices.begin(),
				     actual_dof_indices.end(),
				     new_index);
		  if (*it != new_index)
		    actual_dof_indices.insert(it, new_index);
		}
	    }

	  if (keep_constrained_entries == true)
	    {
	      for (unsigned int j=0; j<n_local_dofs; ++j)
		{
		  sparsity_pattern.add(global_row,
				       local_dof_indices[j]);
		  sparsity_pattern.add(local_dof_indices[j],
				       global_row);
		}
	    }
	  else
	    sparsity_pattern.add(global_row,global_row);
	}

      const unsigned int n_actual_dofs = actual_dof_indices.size();

				   // now add the indices we collected above
				   // to the sparsity pattern. Very easy
				   // here - just add the same array to all
				   // the columns...
      for (unsigned int i=0; i<n_actual_dofs; ++i)
	sparsity_pattern.add_entries(actual_dof_indices[i],
				     actual_dof_indices.begin(),
				     actual_dof_indices.end(),
				     true);
      return;
    }


				   // complicated case: we need to filter
				   // out some indices. then the function
				   // gets similar to the function for
				   // distributing matrix entries, see there
				   // for additional comments.
  std::vector<internals::distributing> my_indices (n_local_dofs);
  std::vector<std::pair<unsigned int, const ConstraintLine *> > constraint_lines;

				   // cache whether we have to resolve any
				   // indirect rows generated from resolving
				   // constrained dofs.
  bool have_indirect_rows = false;
  {
    unsigned int added_rows = 0;
				   // first add the indices in an unsorted
				   // way and only keep track of the
				   // constraints that appear. They are
				   // resolved in a second step.
    for (unsigned int i = 0; i<n_local_dofs; ++i)
      {
	if (constraint_line_exists.size() <= local_dof_indices[i] ||
	    constraint_line_exists[local_dof_indices[i]] == false)
	  {
	    my_indices[added_rows].global_row = local_dof_indices[i];
	    my_indices[added_rows].local_row = i;
	    ++added_rows;
	    continue;
	  }

	ConstraintLine index_comparison;
	index_comparison.line = local_dof_indices[i];

	const std::vector<ConstraintLine>::const_iterator
	  position = std::lower_bound (lines.begin(),
				       lines.end(),
				       index_comparison);
	Assert (position->line == local_dof_indices[i],
		ExcInternalError());

	constraint_lines.push_back (std::make_pair<unsigned int,
				    const ConstraintLine *>(i,&*position));
      }
    Assert (constraint_lines.size() + added_rows == n_local_dofs,
	    ExcInternalError());
    my_indices.resize (added_rows);
  }
  internals::list_shellsort (my_indices);

  Threads::ThreadMutex::ScopedLock lock(mutex);

				   // now in the second step actually
				   // resolve the constraints
  const unsigned int n_constrained_dofs = constraint_lines.size();
  for (unsigned int i=0; i<n_constrained_dofs; ++i)
    {
      const unsigned int local_row = constraint_lines[i].first;
      const unsigned int global_row = local_dof_indices[local_row];
      const ConstraintLine * position = constraint_lines[i].second;
      for (unsigned int q=0; q<position->entries.size(); ++q)
	{
	  have_indirect_rows = true;
	  internals::insert_index(my_indices, position->entries[q].first,
				  std::make_pair<unsigned int,double> 
				  (local_row, position->entries[q].second));
	}

                                   // need to add the whole row and column
                                   // structure in case we keep constrained
                                   // entries. Unfortunately, we can't use
                                   // the nice matrix structure we use
                                   // elsewhere, so manually add those
                                   // indices one by one.
      if (keep_constrained_entries == true)
	{
	  for (unsigned int j=0; j<n_local_dofs; ++j)
	    {
	      if (dof_mask[local_row][j] == true)
		sparsity_pattern.add(global_row,
				     local_dof_indices[j]);
	      if (dof_mask[j][local_row] == true)
		sparsity_pattern.add(local_dof_indices[j],
				     global_row);
	    }
	}
      else
				   // don't keep constrained entries - just
				   // add the diagonal.
	sparsity_pattern.add(global_row,global_row);
    }

  const unsigned int n_actual_dofs = my_indices.size();

				   // create arrays for the column indices
				   // that will then be written into the
				   // sparsity pattern.
  std::vector<unsigned int> cols (n_actual_dofs);

  for (unsigned int i=0; i<n_actual_dofs; ++i)
    {
      std::vector<unsigned int>::iterator col_ptr = cols.begin();
      const unsigned int row = my_indices[i].global_row;
      const unsigned int loc_row = my_indices[i].local_row;

				   // fast function if there are no indirect
				   // references to any of the local rows at
				   // all on this set of dofs
      if (have_indirect_rows == false)
	{
	  Assert(loc_row < n_local_dofs,
		 ExcInternalError());

	  for (unsigned int j=0; j < n_actual_dofs; ++j)
	    {
	      const unsigned int loc_col = my_indices[j].local_row;
	      Assert(loc_col < n_local_dofs, ExcInternalError());

	      if (dof_mask[loc_row][loc_col] == true)
		*col_ptr++ = my_indices[j].global_row;
	    }
	}

				   // slower functions when there are
				   // indirect references and when we need
				   // to do some more checks.
      else
	{
	  for (unsigned int j=0; j < n_actual_dofs; ++j)
	    {
	      const unsigned int loc_col = my_indices[j].local_row;

	      bool add_this = false;

				   // case 1: row has direct contribution in
				   // local matrix
	      if (loc_row != deal_II_numbers::invalid_unsigned_int)
		{
		  Assert (loc_row < n_local_dofs, ExcInternalError());

				   // case 1a: col has direct contribution
				   // in local matrix
		  if (loc_col != deal_II_numbers::invalid_unsigned_int)
		    {
		      Assert (loc_col < n_local_dofs, ExcInternalError());
		      if (dof_mask[loc_row][loc_col] == true)
			goto add_this_index;
		    }

				   // account for indirect contributions by
				   // constraints
		  if (my_indices[j].constraints != 0)
		    {
		      std::vector<std::pair<unsigned int,double> > &constraint_j = 
			*my_indices[j].constraints;

		      for (unsigned int p=0; p<constraint_j.size(); ++p)
			if (dof_mask[loc_row][constraint_j[p].first] == true)
			  goto add_this_index;
		    }
		}

				   // account for indirect contributions by
				   // constraints in row, going trough the
				   // direct and indirect references in the
				   // given column.
	      if (my_indices[i].constraints != 0)
		{
		  std::vector<std::pair<unsigned int,double> > &constraint_i = 
		    *my_indices[i].constraints;
		  for (unsigned int q=0; q<constraint_i.size(); ++q)
		    {
		      if (loc_col != deal_II_numbers::invalid_unsigned_int)
			{
			  Assert (loc_col < n_local_dofs, ExcInternalError());
			  if (dof_mask[constraint_i[q].first][loc_col] == true)
			    goto add_this_index;
			}
  
		      if (my_indices[j].constraints != 0)
			{
			  std::vector<std::pair<unsigned int,double> > &constraint_j = 
			    *my_indices[j].constraints;

			  for (unsigned int p=0; p<constraint_j.size(); ++p)
			    if (dof_mask[constraint_i[q].first]
			                [constraint_j[p].first] == true)
			      goto add_this_index;
			}
		    }
		}

				   // if we got some nontrivial value,
				   // append it to the array of values.
	      if (add_this == true)
		{
		add_this_index:
		  *col_ptr++ = my_indices[j].global_row;
		}
	    }
	}

				   // finally, write all the information
				   // that accumulated under the given
				   // process into the global matrix row and
				   // into the vector
      if (col_ptr != cols.begin())
	sparsity_pattern.add_entries(row, cols.begin(), col_ptr,
				     true);
    }
}




template <typename SparsityType>
inline
void
ConstraintMatrix::
add_entries_local_to_global (const std::vector<unsigned int> &local_dof_indices,
			     SparsityType                    &sparsity_pattern,
			     const bool                       keep_constrained_entries,
			     const Table<2,bool>             &dof_mask,
			     internal::bool2type<true> ) const
{
				   // just as the other
				   // add_entries_local_to_global function,
				   // but now specialized for block
				   // matrices.
  Assert (sparsity_pattern.n_rows() == sparsity_pattern.n_cols(), ExcNotQuadratic());
  Assert (sparsity_pattern.n_block_rows() == sparsity_pattern.n_block_cols(),
	  ExcNotQuadratic());

  const unsigned int n_local_dofs = local_dof_indices.size();
  const unsigned int num_blocks = sparsity_pattern.n_block_rows();

  bool dof_mask_is_active = false;
  if (dof_mask.n_rows() == n_local_dofs)
    {
      dof_mask_is_active = true;
      Assert (dof_mask.n_cols() == n_local_dofs,
	      ExcDimensionMismatch(dof_mask.n_cols(), n_local_dofs));
    }

				   // if the dof mask is not active, all we
				   // have to do is to add some indices in a
				   // matrix format. To do this, we first
				   // create an array of all the indices
				   // that are to be added. these indices
				   // are the local dof indices plus some
				   // indices that come from constraints.
  if (dof_mask_is_active == false)
    {
      std::vector<unsigned int> actual_dof_indices (n_local_dofs);
      unsigned int added_rows = 0;
      bool have_indirect_rows = false;
      std::vector<std::pair<unsigned int, const ConstraintLine *> > constraint_lines;
      for (unsigned int i = 0; i<n_local_dofs; ++i)
	{
	  if (constraint_line_exists.size() <= local_dof_indices[i] ||
	      constraint_line_exists[local_dof_indices[i]] == false)
	    {
	      actual_dof_indices[added_rows] = local_dof_indices[i];
	      ++added_rows;
	      continue;
	    }

	  ConstraintLine index_comparison;
	  index_comparison.line = local_dof_indices[i];

	  const std::vector<ConstraintLine>::const_iterator
	    position = std::lower_bound (lines.begin(),
					 lines.end(),
					 index_comparison);
	  Assert (position->line == local_dof_indices[i],
		  ExcInternalError());

	  constraint_lines.push_back (std::make_pair<unsigned int,
				      const ConstraintLine *>(i,&*position));
      }
      Assert (constraint_lines.size() + added_rows == n_local_dofs,
	      ExcInternalError());
      actual_dof_indices.resize (added_rows);
      std::sort (actual_dof_indices.begin(), actual_dof_indices.end());

      Threads::ThreadMutex::ScopedLock lock(mutex);

      const unsigned int n_constrained_dofs = constraint_lines.size();
      for (unsigned int i=0; i<n_constrained_dofs; ++i)
	{
	  const unsigned int local_row = constraint_lines[i].first;
	  const unsigned int global_row = local_dof_indices[local_row];
	  const ConstraintLine * position = constraint_lines[i].second;
	  for (unsigned int q=0; q<position->entries.size(); ++q)
	    {
	      have_indirect_rows = true;
	      const unsigned int new_index = position->entries[q].first;
	      if (actual_dof_indices.back() < new_index)
		{
		  actual_dof_indices.push_back(new_index);
		}
	      else
		{
		  std::vector<unsigned int>::iterator it = 
		    std::lower_bound(actual_dof_indices.begin(),
				     actual_dof_indices.end(),
				     new_index);
		  if (*it != new_index)
		    actual_dof_indices.insert(it, new_index);
		}
	    }

	  if (keep_constrained_entries == true)
	    {
	      for (unsigned int j=0; j<n_local_dofs; ++j)
		{
		  sparsity_pattern.add(global_row,
				       local_dof_indices[j]);
		  sparsity_pattern.add(local_dof_indices[j],
				       global_row);
		}
	    }
	  else
	    sparsity_pattern.add(global_row,global_row);
	}

      const unsigned int n_actual_dofs = actual_dof_indices.size();

				   // additional construct that also takes
				   // care of block indices.
      std::vector<unsigned int> block_starts(num_blocks+1, n_actual_dofs);
      internals::make_block_starts (sparsity_pattern, actual_dof_indices,
				    block_starts);

				   // easy operation - just go trough the
				   // individual blocks and add the same
				   // array for each row
      for (unsigned int block=0; block<num_blocks; ++block)
	{
	  const unsigned int next_block = block_starts[block+1];
	  for (unsigned int i=block_starts[block]; i<next_block; ++i)
	    {
	      Assert (i<n_actual_dofs, ExcInternalError());
	      const unsigned int row = actual_dof_indices[i];
	      Assert (row < sparsity_pattern.block(block,0).n_rows(),
		      ExcInternalError());
	      std::vector<unsigned int>::iterator index_it = actual_dof_indices.begin();
	      for (unsigned int block_col = 0; block_col<num_blocks; ++block_col)
		{
		  const unsigned int next_block_col = block_starts[block_col+1];
		  sparsity_pattern.block(block,block_col).
		    add_entries(row,
				index_it, 
				actual_dof_indices.begin() + next_block_col, 
				true);
		  index_it = actual_dof_indices.begin() + next_block_col;
		}
	    }
	}
      return;
    }

				   // difficult case with dof_mask, similar
				   // to the distribute_local_to_global
				   // function for block matrices
  std::vector<internals::distributing> my_indices (n_local_dofs);
  std::vector<std::pair<unsigned int, const ConstraintLine *> > constraint_lines;

				   // cache whether we have to resolve any
				   // indirect rows generated from resolving
				   // constrained dofs.
  bool have_indirect_rows = false;
  {
    unsigned int added_rows = 0;
				   // first add the indices in an unsorted
				   // way and only keep track of the
				   // constraints that appear. They are
				   // resolved in a second step.
    for (unsigned int i = 0; i<n_local_dofs; ++i)
      {
	if (constraint_line_exists.size() <= local_dof_indices[i] ||
	    constraint_line_exists[local_dof_indices[i]] == false)
	  {
	    my_indices[added_rows].global_row = local_dof_indices[i];
	    my_indices[added_rows].local_row = i;
	    ++added_rows;
	    continue;
	  }

	ConstraintLine index_comparison;
	index_comparison.line = local_dof_indices[i];

	const std::vector<ConstraintLine>::const_iterator
	  position = std::lower_bound (lines.begin(),
				       lines.end(),
				       index_comparison);
	Assert (position->line == local_dof_indices[i],
		ExcInternalError());

	constraint_lines.push_back (std::make_pair<unsigned int,
				    const ConstraintLine *>(i,&*position));
      }
    Assert (constraint_lines.size() + added_rows == n_local_dofs,
	    ExcInternalError());
    my_indices.resize (added_rows);
  }
  internals::list_shellsort (my_indices);

  Threads::ThreadMutex::ScopedLock lock(mutex);

				   // now in the second step actually
				   // resolve the constraints
  const unsigned int n_constrained_dofs = constraint_lines.size();
  for (unsigned int i=0; i<n_constrained_dofs; ++i)
    {
      const unsigned int local_row = constraint_lines[i].first;
      const unsigned int global_row = local_dof_indices[local_row];
      const ConstraintLine * position = constraint_lines[i].second;
      for (unsigned int q=0; q<position->entries.size(); ++q)
	{
	  have_indirect_rows = true;
	  internals::insert_index(my_indices, position->entries[q].first,
				  std::make_pair<unsigned int,double> 
				  (local_row, position->entries[q].second));
	}

      if (keep_constrained_entries == true)
	{
	  for (unsigned int j=0; j<n_local_dofs; ++j)
	    {
	      if (dof_mask[local_row][j] == true)
		sparsity_pattern.add(global_row,
				     local_dof_indices[j]);
	      if (dof_mask[j][local_row] == true)
		sparsity_pattern.add(local_dof_indices[j],
				     global_row);
	    }
	}
      else
				   // don't keep constrained entries - just
				   // add the diagonal.
	sparsity_pattern.add(global_row,global_row);
    }


  const unsigned int n_actual_dofs = my_indices.size();

  std::vector<unsigned int> localized_indices (n_actual_dofs);
  for (unsigned int i=0; i<n_actual_dofs; ++i)
    localized_indices[i] = my_indices[i].global_row;

				   // additional construct that also takes
				   // care of block indices.
  std::vector<unsigned int> block_starts(num_blocks+1, n_actual_dofs);
  internals::make_block_starts(sparsity_pattern, localized_indices,
			       block_starts);

  std::vector<unsigned int> cols (n_actual_dofs);

				   // the basic difference to the
				   // non-block variant from now onwards
				   // is that we go through the blocks
				   // of the matrix separately.
  for (unsigned int block=0; block<num_blocks; ++block)
    {
      const unsigned int next_block = block_starts[block+1];
      for (unsigned int i=block_starts[block]; i<next_block; ++i)
	{
	  const unsigned int row = localized_indices[i];
	  const unsigned int loc_row = my_indices[i].local_row;

	  for (unsigned int block_col=0; block_col<num_blocks; ++block_col)
	    {
	      const unsigned int next_block_col = block_starts[block_col+1];
	      std::vector<unsigned int>::iterator col_ptr = cols.begin();
	      if (have_indirect_rows == false)
		{
		  Assert(loc_row < n_local_dofs,
			 ExcInternalError());

		  for (unsigned int j=block_starts[block_col]; j < next_block_col; ++j)
		    {
		      const unsigned int loc_col = my_indices[j].local_row;
		      Assert(loc_col < n_local_dofs,
			     ExcInternalError());

		      if (dof_mask[loc_row][loc_col] == true)
			*col_ptr++ = localized_indices[j];
		    }
		}

				   // have indirect references by
				   // constraints, resolve them
	      else
		{
		  for (unsigned int j=block_starts[block_col]; j < next_block_col; ++j)
		    {
		      const unsigned int loc_col = my_indices[j].local_row;

		      bool add_this = false;

		      if (loc_row != deal_II_numbers::invalid_unsigned_int)
			{
			  Assert (loc_row < n_local_dofs,
				  ExcInternalError());

			  if (loc_col != deal_II_numbers::invalid_unsigned_int)
			    {
			      Assert (loc_col < n_local_dofs,
				      ExcInternalError());
			      if (dof_mask[loc_row][loc_col] == true)
				goto add_this_index;
			    }

				   // account for indirect contributions by
				   // constraints
			  if (my_indices[j].constraints != 0)
			    {
			      std::vector<std::pair<unsigned int,double> > 
				&constraint_j = *my_indices[j].constraints;

			      for (unsigned int p=0; p<constraint_j.size(); ++p)
				if (dof_mask[loc_row][constraint_j[p].first] == true)
				  goto add_this_index;
			    }
			}

				   // account for indirect contributions by
				   // constraints in row, going trough the
				   // direct and indirect references in the
				   // given column.
		      if (my_indices[i].constraints != 0)
			{
			  std::vector<std::pair<unsigned int,double> > 
			    &constraint_i = *my_indices[i].constraints;
			  for (unsigned int q=0; q<constraint_i.size(); ++q)
			    {
			      if (loc_col != deal_II_numbers::invalid_unsigned_int)
				{
				  Assert (loc_col < n_local_dofs,
					  ExcInternalError());
				  if (dof_mask[constraint_i[q].first][loc_col] == true)
				    goto add_this_index;
				}
  
			      if (my_indices[j].constraints != 0)
				{
				  std::vector<std::pair<unsigned int,double> > 
				    &constraint_j = *my_indices[j].constraints;

				  for (unsigned int p=0; p<constraint_j.size(); ++p)
				    if (dof_mask[constraint_i[q].first]
			                        [constraint_j[p].first] == true)
				      goto add_this_index;
				}
			    }
			}
		      if (add_this == true)
			{
			add_this_index:
			  *col_ptr++ = localized_indices[j];
			}
		    }
		}

				   // finally, write all the information
				   // that accumulated under the given
				   // process into the global matrix row and
				   // into the vector
	      sparsity_pattern.block(block, block_col).add_entries(row, 
								   cols.begin(), 
								   col_ptr,
								   true);
	    }
	}
    }
}


DEAL_II_NAMESPACE_CLOSE

#endif
