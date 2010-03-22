//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__constraint_matrix_templates_h
#define __deal2__constraint_matrix_templates_h


#include <lac/constraint_matrix.h>

#include <base/table.h>
#include <lac/full_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
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
          += (static_cast<typename VectorType::value_type>
	      (vec(constraint_line->line)) *
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
                            VectorType                      &global_vector,
			    const FullMatrix<double>        &local_matrix) const
{
  Assert (local_vector.size() == local_dof_indices.size(),
          ExcDimensionMismatch(local_vector.size(), local_dof_indices.size()));
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (local_matrix.m() == local_dof_indices.size(),
	  ExcDimensionMismatch(local_matrix.m(), local_dof_indices.size()));
  Assert (local_matrix.n() == local_dof_indices.size(),
	  ExcDimensionMismatch(local_matrix.n(), local_dof_indices.size()));

  const unsigned int n_local_dofs = local_vector.size();
  if (lines.size() == 0)
    global_vector.add(local_dof_indices, local_vector);
  else
    for (unsigned int i=0; i<n_local_dofs; ++i)
      {
				// check whether the current index is
				// constrained. if not, just write the entry
				// into the vector. otherwise, need to resolve
				// the constraint
	if (is_constrained(local_dof_indices[i]) == false)
	  {
	    global_vector(local_dof_indices[i]) += local_vector(i);
	    continue;
	  }

				// find the constraint line to the given
				// global dof index
	const unsigned int line_index = calculate_line_index (local_dof_indices[i]);
	const ConstraintLine * position =
	  lines_cache.size() <= line_index ? 0 : &lines[lines_cache[line_index]];

				// Gauss elimination of the matrix columns
				// with the inhomogeneity. Go through them one
				// by one and again check whether they are
				// constrained. If so, distribute the constraint
	const double val = position->inhomogeneity;
	if (val != 0)
	  for (unsigned int j=0; j<n_local_dofs; ++j)
	    if (is_constrained(local_dof_indices[j]) == false)
	      global_vector(local_dof_indices[j]) -= val * local_matrix(j,i);
	    else
	      {
		const double matrix_entry = local_matrix(j,i);
		if (matrix_entry == 0)
		  continue;

		const ConstraintLine & position_j =
		  lines[lines_cache[calculate_line_index(local_dof_indices[j])]];
		for (unsigned int q=0; q<position_j.entries.size(); ++q)
		  {
		    Assert (is_constrained(position_j.entries[q].first) == false,
			    ExcMessage ("Tried to distribute to a fixed dof."));
		    global_vector(position_j.entries[q].first)
		      -= val * position_j.entries[q].second * matrix_entry;
		  }
	      }

				   // now distribute the constraint,
				   // but make sure we don't touch
				   // the entries of fixed dofs
	for (unsigned int j=0; j<position->entries.size(); ++j)
	  {
	    Assert (is_constrained(position->entries[j].first) == false,
		    ExcMessage ("Tried to distribute to a fixed dof."));
	    global_vector(position->entries[j].first)
	      += local_vector(i) * position->entries[j].second;
	  }
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
      typename VectorType::value_type
	new_value = next_constraint->inhomogeneity;
      for (unsigned int i=0; i<next_constraint->entries.size(); ++i)
	new_value += (static_cast<typename VectorType::value_type>
		      (vec(next_constraint->entries[i].first)) *
                      next_constraint->entries[i].second);
      vec(next_constraint->line) = new_value;
    }
}



				   // Some helper definitions for the
				   // local_to_global functions.
namespace internals
{
				   // this struct contains all the information
				   // we need to store about each of the
				   // global entries (global_row): are they
				   // obtained directly by some local entry
				   // (local_row) or some constraints
				   // (constraint_position). This is not
				   // directly used in
				   // the user code, but accessed via the
				   // GlobalRowsFromLocal.
  struct Distributing
  {
    Distributing (const unsigned int global_row = numbers::invalid_unsigned_int,
		  const unsigned int local_row = numbers::invalid_unsigned_int);
    Distributing (const Distributing &in);
    Distributing & operator = (const Distributing &in);
    bool operator < (const Distributing &in) const {return global_row<in.global_row;};

    unsigned int global_row;
    unsigned int local_row;
    mutable unsigned int constraint_position;
  };

  inline
  Distributing::Distributing (const unsigned int global_row,
			      const unsigned int local_row) :
    global_row (global_row),
    local_row (local_row),
    constraint_position (numbers::invalid_unsigned_int) {}

  inline
  Distributing::Distributing (const Distributing &in) :
    constraint_position (numbers::invalid_unsigned_int)
  {*this = (in);}

  inline
  Distributing & Distributing::operator = (const Distributing &in)
  {
    global_row = in.global_row;
    local_row = in.local_row;
				   // the constraints pointer should not
				   // contain any data here.
    Assert (constraint_position == numbers::invalid_unsigned_int, 
	    ExcInternalError());

    if (in.constraint_position != numbers::invalid_unsigned_int)
      {
	constraint_position = in.constraint_position;
	in.constraint_position = numbers::invalid_unsigned_int;
      }
    return *this;
  }



				    // this is a cache for constraints that
				    // are encountered on a local level.
				    // corresponds to functionality also
				    // provided by
				    // std::vector<std::vector<std::pair<uint,double>
				    // > >, but tuned so that frequent memory
				    // allocation for each entry is
				    // avoided. This is not directly used in
				    // the user code, but accessed via the
				    // GlobalRowsFromLocal.
  struct DataCache
  {
      DataCache ()
		      :
		      element_size (0),
		      data (0),
		      n_used_elements(0)
	{}

      ~DataCache()
	{
	  if (data != 0)
	    delete [] data;
	}

      void reinit ()
	{
	  Assert (element_size == 0, ExcInternalError());
	  element_size = 6;
	  data = new std::pair<unsigned int,double> [20*6];
	  individual_size.resize(20);
	  n_used_elements = 0;
	}

      unsigned int element_size;

      std::pair<unsigned int,double> * data;

      std::vector<unsigned int> individual_size;

      unsigned int n_used_elements;

      unsigned int insert_new_index (const std::pair<unsigned int,double> &pair)
	{
	  if (element_size == 0)
	    reinit();
	  if (n_used_elements == individual_size.size())
	    {
	      std::pair<unsigned int,double> * new_data =
		new std::pair<unsigned int,double> [2*individual_size.size()*element_size];
	      memcpy (new_data, data, individual_size.size()*element_size*
		      sizeof(std::pair<unsigned int,double>));
	      delete [] data;
	      data = new_data;
	      individual_size.resize (2*individual_size.size(), 0);
	    }
	  unsigned int index = n_used_elements;
	  data[index*element_size] = pair;
	  individual_size[index] = 1;
	  ++n_used_elements;
	  return index;
	}

      void append_index (const unsigned int index,
			 const std::pair<unsigned int,double> &pair)
	{
	  Assert (index < n_used_elements, ExcIndexRange (index, 0, n_used_elements));
	  const unsigned int my_size = individual_size[index];
	  if (my_size == element_size)
	    {
	      std::pair<unsigned int,double> * new_data =
		new std::pair<unsigned int,double> [2*individual_size.size()*element_size];
	      for (unsigned int i=0; i<n_used_elements; ++i)
		memcpy (&new_data[i*element_size*2], &data[i*element_size],
			element_size*sizeof(std::pair<unsigned int,double>));
	      delete [] data;
	      data = new_data;
	      element_size *= 2;
	    }
	  data[index*element_size+my_size] = pair;
	  individual_size[index]++;
	}

      unsigned int
      get_size (const unsigned int index) const
	{
	  return individual_size[index];
	}

      const std::pair<unsigned int,double> *
      get_entry (const unsigned int index) const
	{
	  return &data[index*element_size];
	}
  };



				// collects all the global rows from a
				// local contribution (cell) and their
				// origin (direct/constraint). this is
				// basically a vector consisting of
				// "Distributing" structs using access via
				// the DataCache. Provides some
				// specialized sort and insert functions.
				// 
				// in case there are no constraints, this is
				// basically a list of pairs <uint,unit> with
				// the first index being the global index and
				// the second index the local index. The list
				// is sorted with respect to the global index.
				// 
				// in case there are constraints, a global dof
				// might get a contribution also because it
				// gets data from a constrained dof. This
				// means that a global dof might also have
				// indirect contributions from a local dof via
				// a constraint, besides the direct ones.
  struct GlobalRowsFromLocal
  {
    GlobalRowsFromLocal (const unsigned int n_local_rows)
      :
      total_row_indices (n_local_rows), n_active_rows (n_local_rows),
      n_inhomogeneous_rows (0){};

				// implemented below
    void insert_index (const unsigned int global_row,
		       const unsigned int local_row,
		       const double       constraint_value);
    void sort ();

				// return all kind of information on the
				// constraints

				// returns the number of global indices in the
				// struct
    unsigned int size () const { return n_active_rows; };

				// returns the global index of the
				// counter_index-th entry in the list
    unsigned int & global_row (const unsigned int counter_index)
      { return total_row_indices[counter_index].global_row; };

				// returns the number of constraints that are
				// associated to the counter_index-th entry in
				// the list
    unsigned int size (const unsigned int counter_index) const
      { return (total_row_indices[counter_index].constraint_position ==
		numbers::invalid_unsigned_int ?
		0 :
		data_cache.get_size(total_row_indices[counter_index].
				    constraint_position)); };

				// returns the global row associated with the
				// counter_index-th entry in the list
    const unsigned int & global_row (const unsigned int counter_index) const
      { return total_row_indices[counter_index].global_row; };

				// returns the local row in the cell matrix
				// associated with the counter_index-th entry
				// in the list. Returns invalid_unsigned_int
				// for invalid unsigned ints
    const unsigned int & local_row (const unsigned int counter_index) const
      { return total_row_indices[counter_index].local_row; };

				// writable index
    unsigned int & local_row (const unsigned int counter_index)
      { return total_row_indices[counter_index].local_row; };

				// returns the local row in the cell matrix
				// associated with the counter_index-th entry
				// in the list in the index_in_constraint-th
				// position of constraints
    unsigned int local_row (const unsigned int counter_index,
			    const unsigned int index_in_constraint) const
      { return (data_cache.get_entry(total_row_indices[counter_index].constraint_position)
		[index_in_constraint]).first; };

				// returns the value of the constraint in the
				// counter_index-th entry in the list in the
				// index_in_constraint-th position of
				// constraints
    double constraint_value (const unsigned int counter_index,
			     const unsigned int index_in_constraint) const
      { return (data_cache.get_entry(total_row_indices[counter_index].constraint_position)
		[index_in_constraint]).second; };

				// returns whether there is one row with
				// indirect contributions (i.e., there has
				// been at least one constraint with
				// non-trivial ConstraintLine)
    bool have_indirect_rows () const { return data_cache.element_size; };

				// append an entry that is constrained. This
				// means that there is one less row that data
				// needs to be inserted into.
    void insert_constraint (const unsigned int constrained_local_dof)
      { --n_active_rows;
	total_row_indices[n_active_rows].local_row = constrained_local_dof; }

				// last row that was set to be constrained
    unsigned int last_constrained_local_row ()
      { Assert (total_row_indices.back().global_row == numbers::invalid_unsigned_int,
		ExcInternalError());
	return total_row_indices.back().local_row; };

				// remove last entry in the list of global
				// rows. Some elements at the end of the list
				// total_row_indices are used to temporarily
				// collect constrained dofs, but they should
				// not have a global row index. Check that and
				// eliminate the last such entry.
    void pop_back ()
      { Assert (total_row_indices.back().global_row == numbers::invalid_unsigned_int,
		ExcInternalError());
	total_row_indices.pop_back(); };

				// store that the constraint_dof at the last
				// position is inhomogeneous and needs to be
				// considered further.
    void set_last_row_inhomogeneous ()
      { std::swap (total_row_indices[n_active_rows+n_inhomogeneous_rows].local_row,
		   total_row_indices.back().local_row);
	++n_inhomogeneous_rows; };

				// 
    unsigned int inhomogeneity (unsigned int i) const
      { return total_row_indices[n_active_rows+i].local_row; };

				// a vector that contains all the global ids
				// and the corresponding local ids as well as
				// a pointer to that data where we store how
				// to resolve constraints.
    std::vector<Distributing> total_row_indices;

				// holds the actual data from the constraints
    DataCache                 data_cache;

				// how many rows there are
    unsigned int              n_active_rows;

				// the number of rows with inhomogeneous
				// constraints
    unsigned int              n_inhomogeneous_rows;
  };

				   // a function that appends an additional
				   // row to the list of values, or appends a
				   // value to an already existing
				   // row. Similar functionality as for
				   // std::map<unsigned int,Distributing>, but
				   // here done for a
				   // std::vector<Distributing>, much faster
				   // for short lists as we have them here
  inline
  void
  GlobalRowsFromLocal::insert_index (const unsigned int global_row,
				     const unsigned int local_row,
				     const double       constraint_value)
  {
    typedef std::vector<Distributing>::iterator index_iterator;
    index_iterator pos, pos1;
    Distributing row_value (global_row);
    std::pair<unsigned int,double> constraint (local_row, constraint_value);

				   // check whether the list was really
				   // sorted before entering here
    for (unsigned int i=1; i<n_active_rows; ++i)
      Assert (total_row_indices[i-1] < total_row_indices[i], ExcInternalError());

    pos = std::lower_bound (total_row_indices.begin(),
			    total_row_indices.begin()+n_active_rows,
			    row_value);
    if (pos->global_row == global_row)
      pos1 = pos;
    else
      {
	pos1 = total_row_indices.insert(pos, row_value);
	++n_active_rows;
      }

    if (pos1->constraint_position == numbers::invalid_unsigned_int)
      pos1->constraint_position = data_cache.insert_new_index (constraint);
    else
      data_cache.append_index (pos1->constraint_position, constraint);
  }

				   // this sort algorithm sorts
				   // std::vector<Distributing>, but does not
				   // take the constraints into account. this
				   // means that in case that constraints are
				   // already inserted, this function does not
				   // work as expected. Use shellsort, which
				   // is very fast in case the indices are
				   // already sorted (which is the usual case
				   // with DG elements), and not too slow in
				   // other cases
  inline
  void
  GlobalRowsFromLocal::sort ()
  {
    unsigned int i, j, j2, temp, templ, istep;
    unsigned int step;

				   // check whether the
				   // constraints are really empty.
    const unsigned int length = size();
#ifdef DEBUG
    for (unsigned int i=0; i<length; ++i)
      Assert (total_row_indices[i].constraint_position ==
	      numbers::invalid_unsigned_int,
	      ExcInternalError());
#endif

    step = length/2;
    while (step > 0)
      {
	for (i=step; i < length; i++)
	  {
	    istep = step;
	    j = i;
	    j2 = j-istep;
	    temp = total_row_indices[i].global_row;
	    templ = total_row_indices[i].local_row;
	    if (total_row_indices[j2].global_row > temp)
	      {
		while ((j >= istep) && (total_row_indices[j2].global_row > temp))
		  {
		    total_row_indices[j].global_row = total_row_indices[j2].global_row;
		    total_row_indices[j].local_row = total_row_indices[j2].local_row;
		    j = j2;
		    j2 -= istep;
		  }
		total_row_indices[j].global_row = temp;
		total_row_indices[j].local_row = templ;
	      }
	  }
	step = step>>1;
      }
  }

				   // function for block matrices: Find out
				   // where in the list of local dofs (sorted
				   // according to global ids) the individual
				   // blocks start. Transform the global
				   // indices to block-local indices in order
				   // to be able to use functions like
				   // vector.block(1)(block_local_id), instead
				   // of vector(global_id). This avoids
				   // transforming indices one-by-one later
				   // on.
  template <class BlockType>
  inline
  void
  make_block_starts (const BlockType           &block_object,
		     GlobalRowsFromLocal       &global_rows,
		     std::vector<unsigned int> &block_starts)
  {
    Assert (block_starts.size() == block_object.n_block_rows() + 1,
	    ExcDimensionMismatch(block_starts.size(),
				 block_object.n_block_rows()+1));

    typedef std::vector<Distributing>::iterator row_iterator;
    row_iterator block_indices = global_rows.total_row_indices.begin();

    const unsigned int num_blocks = block_object.n_block_rows();
    const unsigned int n_active_rows = global_rows.size();

				   // find end of rows.
    block_starts[0] = 0;
    for (unsigned int i=1;i<num_blocks;++i)
      {
	row_iterator first_block =
	  std::lower_bound (block_indices,
			    global_rows.total_row_indices.begin()+n_active_rows,
			    Distributing(block_object.get_row_indices().block_start(i)));
	block_starts[i] = first_block - global_rows.total_row_indices.begin();
	block_indices = first_block;
      }

				   // transform row indices to block-local
				   // index space
    for (unsigned int i=block_starts[1]; i<n_active_rows; ++i)
      global_rows.global_row(i) = block_object.get_row_indices().
	global_to_local(global_rows.global_row(i)).second;
  }



				   // same as before, but for
				   // std::vector<uint> instead of
				   // GlobalRowsFromLocal. Used in functions
				   // for sparsity patterns.
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
      row_indices[i] = block_object.get_row_indices().
	global_to_local(row_indices[i]).second;
  }



				   // resolves constraints of one column at
				   // the innermost loop. goes through the
				   // origin of each global entry and finds
				   // out which data we need to collect.
  inline
  double resolve_matrix_entry (const GlobalRowsFromLocal&global_rows,
			       const unsigned int        i,
			       const unsigned int        j,
			       const unsigned int        loc_row,
			       const FullMatrix<double> &local_matrix)
  {
    const unsigned int loc_col = global_rows.local_row(j);
    double col_val;

				   // case 1: row has direct contribution in
				   // local matrix. decide whether col has a
				   // direct contribution. if not,
				   // set the value to zero.
    if (loc_row != numbers::invalid_unsigned_int)
      {
	col_val = ((loc_col != numbers::invalid_unsigned_int) ? 
		   local_matrix(loc_row, loc_col) : 0);

				   // account for indirect contributions by
				   // constraints in column
	for (unsigned int p=0; p<global_rows.size(j); ++p)
	  col_val += (local_matrix(loc_row,global_rows.local_row(j,p)) *
		      global_rows.constraint_value(j,p));
      }

				   // case 2: row has no direct contribution in
				   // local matrix
    else
      col_val = 0;

				   // account for indirect contributions by
				   // constraints in row, going trough the
				   // direct and indirect references in the
				   // given column.
    for (unsigned int q=0; q<global_rows.size(i); ++q)
      {
	double add_this = loc_col != numbers::invalid_unsigned_int ?
	  local_matrix(global_rows.local_row(i,q), loc_col) : 0;

	for (unsigned int p=0; p<global_rows.size(j); ++p)
	  add_this += (local_matrix(global_rows.local_row(i,q),
				    global_rows.local_row(j,p))
		       *
		       global_rows.constraint_value(j,p));
	col_val += add_this * global_rows.constraint_value(i,q);
      }
    return col_val;
  }



				// computes all entries that need to be
				// written into global_rows[i]. Lists the
				// resulting values in val_ptr, and the
				// corresponding column indices in col_ptr.
  template <typename number>
  inline
  void
  resolve_matrix_row (const GlobalRowsFromLocal&global_rows,
		      const unsigned int        i,
		      const unsigned int        column_start,
		      const unsigned int        column_end,
		      const FullMatrix<double> &local_matrix,
		      unsigned int *           &col_ptr,
		      number *                 &val_ptr)
  {
    Assert (global_rows.size() >= column_end,
	    ExcIndexRange (column_end, 0, global_rows.size()));
    const unsigned int loc_row = global_rows.local_row(i);

				   // fast function if there are no indirect
				   // references to any of the local rows at
				   // all on this set of dofs (saves a lot
				   // of checks). the only check we actually
				   // need to perform is whether the matrix
				   // element is zero.
    if (global_rows.have_indirect_rows() == false)
      {
	Assert(loc_row < local_matrix.m(), ExcInternalError());
	const double * matrix_ptr = &local_matrix(loc_row, 0);

	for (unsigned int j=column_start; j<column_end; ++j)
	  {
	    const unsigned int loc_col = global_rows.local_row(j);
	    const double col_val = matrix_ptr[loc_col];
	    if (col_val != 0.)
	      {
		*val_ptr++ = static_cast<number> (col_val);
		*col_ptr++ = global_rows.global_row(j);
	      }
	  }
      }

				   // more difficult part when there are
				   // indirect references and when we need
				   // to do some more checks.
    else
      {
	for (unsigned int j=column_start; j<column_end; ++j)
	  {
	    double col_val = resolve_matrix_entry (global_rows, i, j,
						   loc_row, local_matrix);

				   // if we got some nontrivial value,
				   // append it to the array of values.
	    if (col_val != 0.)
	      {
		*val_ptr++ = static_cast<number> (col_val);
		*col_ptr++ = global_rows.global_row(j);
	      }
	  }
      }
  }



				   // specialized function that can write into
				   // the row of a SparseMatrix<number>.
  namespace dealiiSparseMatrix
  {
    template <typename number>
    inline
    void add_value (const double value,
		    const unsigned int row,
		    const unsigned int column,
		    const unsigned int * col_ptr,
		    const bool   are_on_diagonal,
		    unsigned int &counter,
		    number       *val_ptr)
    {
      if (value != 0.)
	{
	  Assert (col_ptr != 0,
		  typename SparseMatrix<number>::ExcInvalidIndex (row, column));
	  if (are_on_diagonal)
	    {
	      val_ptr[0] += value;
	      return;
	    }
	  while (col_ptr[counter] < column)
	    ++counter;
	  Assert (col_ptr[counter] == column,
		  typename SparseMatrix<number>::ExcInvalidIndex(row, column));
	  val_ptr[counter] += static_cast<number>(value);
	}
    }
  }


				// similar as before, now with shortcut for
				// deal.II sparse matrices. this lets use
				// avoid using extra arrays, and does all the
				// operations just in place, i.e., in the
				// respective matrix row
  template <typename number>
  inline
  void
  resolve_matrix_row (const GlobalRowsFromLocal&global_rows,
		      const unsigned int        i,
		      const unsigned int        column_start,
		      const unsigned int        column_end,
		      const FullMatrix<double> &local_matrix,
		      SparseMatrix<number>     *sparse_matrix)
  {
    Assert (global_rows.size() >= column_end,
	    ExcIndexRange (column_end, 0, global_rows.size()));
    const SparsityPattern & sparsity = sparse_matrix->get_sparsity_pattern();
#ifndef DEBUG
    if (sparsity.n_nonzero_elements() == 0)
      return;
#endif
    const std::size_t * row_start = sparsity.get_rowstart_indices();
    const unsigned int * sparsity_struct = sparsity.get_column_numbers();

    const unsigned int row = global_rows.global_row(i);
    const unsigned int loc_row = global_rows.local_row(i);

#ifdef DEBUG
    const unsigned int * col_ptr = sparsity.n_nonzero_elements() == 0 ? 0 :
      &sparsity_struct[row_start[row]];
    number * val_ptr = sparsity.n_nonzero_elements() == 0 ? 0 :
      &sparse_matrix->global_entry (row_start[row]);
#else
    const unsigned int * col_ptr = &sparsity_struct[row_start[row]];
    number * val_ptr = &sparse_matrix->global_entry (row_start[row]);
#endif
    const bool optimize_diagonal = sparsity.optimize_diagonal();
    unsigned int counter = optimize_diagonal;

				// distinguish three cases about what can
				// happen for checking whether the diagonal is
				// the first element of the row. this avoids
				// if statements at the innermost loop
				// positions

    if (!optimize_diagonal) // case 1: no diagonal optimization in matrix
      {
	if (global_rows.have_indirect_rows() == false)
	  {
	    Assert(loc_row < local_matrix.m(),
		   ExcIndexRange(loc_row, 0, local_matrix.m()));
	    const double * matrix_ptr = &local_matrix(loc_row, 0);

	    for (unsigned int j=column_start; j<column_end; ++j)
	      {
		const unsigned int loc_col = global_rows.local_row(j);
		const double col_val = matrix_ptr[loc_col];
		dealiiSparseMatrix::add_value (col_val, row, 
					       global_rows.global_row(j), 
					       col_ptr, false, counter, 
					       val_ptr);
	      }
	  }
	else
	  {
	    for (unsigned int j=column_start; j<column_end; ++j)
	      {
		double col_val = resolve_matrix_entry (global_rows, i, j,
						       loc_row, local_matrix);
		dealiiSparseMatrix::add_value (col_val, row, 
					       global_rows.global_row(j), col_ptr,
					       false, counter, val_ptr);
	      }
	  }
      }
    else if (i>=column_start && i<column_end) // case 2: can split loop
      {
	if (global_rows.have_indirect_rows() == false)
	  {
	    Assert(loc_row < local_matrix.m(),
		   ExcIndexRange(loc_row, 0, local_matrix.m()));
	    const double * matrix_ptr = &local_matrix(loc_row, 0);

	    for (unsigned int j=column_start; j<i; ++j)
	      {
		const unsigned int loc_col = global_rows.local_row(j);
		const double col_val = matrix_ptr[loc_col];
		dealiiSparseMatrix::add_value(col_val, row, 
					      global_rows.global_row(j), col_ptr,
					      false, counter, val_ptr);
	      }
	    val_ptr[0] += matrix_ptr[loc_row];
	    for (unsigned int j=i+1; j<column_end; ++j)
	      {
		const unsigned int loc_col = global_rows.local_row(j);
		const double col_val = matrix_ptr[loc_col];
		dealiiSparseMatrix::add_value(col_val, row, 
					      global_rows.global_row(j), col_ptr,
					      false, counter, val_ptr);
	      }
	  }
	else
	  {
	    for (unsigned int j=column_start; j<i; ++j)
	      {
		double col_val = resolve_matrix_entry (global_rows, i, j,
						       loc_row, local_matrix);
		dealiiSparseMatrix::add_value (col_val, row, 
					       global_rows.global_row(j), col_ptr,
					       false, counter, val_ptr);
	      }
	    val_ptr[0] += resolve_matrix_entry (global_rows, i, i, loc_row,
						local_matrix);
	    for (unsigned int j=i+1; j<column_end; ++j)
	      {
		double col_val = resolve_matrix_entry (global_rows, i, j,
						       loc_row, local_matrix);
		dealiiSparseMatrix::add_value (col_val, row, 
					       global_rows.global_row(j), col_ptr,
					       false, counter, val_ptr);
	      }
	  }
      }
				// case 3: can't say - need to check inside
				// the loop
    else if (global_rows.have_indirect_rows() == false)
      {
	Assert(loc_row < local_matrix.m(),
	       ExcIndexRange(loc_row, 0, local_matrix.m()));
	const double * matrix_ptr = &local_matrix(loc_row, 0);

	for (unsigned int j=column_start; j<column_end; ++j)
	  {
	    const unsigned int loc_col = global_rows.local_row(j);
	    const double col_val = matrix_ptr[loc_col];
	    dealiiSparseMatrix::add_value(col_val, row, 
					  global_rows.global_row(j), col_ptr,
					  row==global_rows.global_row(j), 
					  counter, val_ptr);
	  }
      }
    else
      {
	for (unsigned int j=column_start; j<column_end; ++j)
	  {
	    double col_val = resolve_matrix_entry (global_rows, i, j,
						   loc_row, local_matrix);
	    dealiiSparseMatrix::add_value (col_val, row, 
					   global_rows.global_row(j), col_ptr,
					   row==global_rows.global_row(j), 
					   counter, val_ptr);
	  }
      }
  }



				// Same function to resolve all entries that
				// will be added to the given global row
				// global_rows[i] as before, now for sparsity
				// pattern
  inline
  void
  resolve_matrix_row (const GlobalRowsFromLocal     &global_rows,
		      const unsigned int             i,
		      const unsigned int             column_start,
		      const unsigned int             column_end,
		      const Table<2,bool>           &dof_mask,
		      std::vector<unsigned int>::iterator &col_ptr)
  {
    const unsigned int loc_row = global_rows.local_row(i);

				   // fast function if there are no indirect
				   // references to any of the local rows at
				   // all on this set of dofs
    if (global_rows.have_indirect_rows() == false)
      {
	Assert(loc_row < dof_mask.n_rows(),
	       ExcInternalError());

	for (unsigned int j=column_start; j<column_end; ++j)
	  {
	    const unsigned int loc_col = global_rows.local_row(j);
	    Assert(loc_col < dof_mask.n_cols(), ExcInternalError());

	    if (dof_mask[loc_row][loc_col] == true)
	      *col_ptr++ = global_rows.global_row(j);
	  }
      }

				   // slower functions when there are
				   // indirect references and when we need
				   // to do some more checks.
    else
      {
	for (unsigned int j=column_start; j<column_end; ++j)
	  {
	    const unsigned int loc_col = global_rows.local_row(j);
	    if (loc_row != numbers::invalid_unsigned_int)
	      {
		Assert (loc_row < dof_mask.n_rows(), ExcInternalError());
		if (loc_col != numbers::invalid_unsigned_int)
		  {
		    Assert (loc_col < dof_mask.n_cols(), ExcInternalError());
		    if (dof_mask[loc_row][loc_col] == true)
		      goto add_this_index;
		  }

		for (unsigned int p=0; p<global_rows.size(j); ++p)
		  if (dof_mask[loc_row][global_rows.local_row(j,p)] == true)
		    goto add_this_index;
	      }

	    for (unsigned int q=0; q<global_rows.size(i); ++q)
	      {
		if (loc_col != numbers::invalid_unsigned_int)
		  {
		    Assert (loc_col < dof_mask.n_cols(), ExcInternalError());
		    if (dof_mask[global_rows.local_row(i,q)][loc_col] == true)
		      goto add_this_index;
		  }

		for (unsigned int p=0; p<global_rows.size(j); ++p)
		  if (dof_mask[global_rows.local_row(i,q)]
		              [global_rows.local_row(j,p)] == true)
		    goto add_this_index;
	      }

	    continue;
				   // if we got some nontrivial value,
				   // append it to the array of values.
	  add_this_index:
	    *col_ptr++ = global_rows.global_row(j);
	  }
      }
  }


} // end of namespace internals



void
ConstraintMatrix::
make_sorted_row_list (const std::vector<unsigned int> &local_dof_indices,
		      internals::GlobalRowsFromLocal  &global_rows) const
{
  const unsigned int n_local_dofs = local_dof_indices.size();
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

				   // cache whether we have to resolve any
				   // indirect rows generated from resolving
				   // constrained dofs.
  unsigned int added_rows = 0;

				   // first add the indices in an unsorted
				   // way and only keep track of the
				   // constraints that appear. They are
				   // resolved in a second step.
  for (unsigned int i = 0; i<n_local_dofs; ++i)
    {
      if (is_constrained(local_dof_indices[i]) == false)
	{
	  global_rows.global_row(added_rows)  = local_dof_indices[i];
	  global_rows.local_row(added_rows++) = i;
	  continue;
	}
      global_rows.insert_constraint(i);
    }
  global_rows.sort();

  const unsigned int n_constrained_rows = n_local_dofs-added_rows;
  for (unsigned int i=0; i<n_constrained_rows; ++i)
    {
      const unsigned int local_row = global_rows.last_constrained_local_row();
      Assert (local_row < n_local_dofs, ExcIndexRange(local_row, 0, n_local_dofs));
      const unsigned int global_row = local_dof_indices[local_row];
      Assert (is_constrained(global_row), ExcInternalError());
      const ConstraintLine & position =
	lines[lines_cache[calculate_line_index(global_row)]];
      if (position.inhomogeneity != 0)
	global_rows.set_last_row_inhomogeneous();
      else
	global_rows.pop_back();
      for (unsigned int q=0; q<position.entries.size(); ++q)
	global_rows.insert_index (position.entries[q].first,
				  local_row,
				  position.entries[q].second);
    }
}




//TODO: This function is DANGEROUS, because it does more than it claims!


				// Basic idea of setting up a list of
				// all global dofs: first find all rows and columns
				// that we are going to write touch,
				// and then go through the
				// lines and collect all the local rows that
				// are related to it.
template <typename MatrixType>
inline
void
ConstraintMatrix::
make_sorted_row_list (const FullMatrix<double>        &local_matrix,
		      const std::vector<unsigned int> &local_dof_indices,
		      MatrixType                      &global_matrix,
		      internals::GlobalRowsFromLocal  &global_rows) const
{
  const unsigned int n_local_dofs = local_dof_indices.size();

  double average_diagonal = 0;
  for (unsigned int i=0; i<n_local_dofs; ++i)
    average_diagonal += std::fabs (local_matrix(i,i));
  average_diagonal /= static_cast<double>(n_local_dofs);

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

				   // cache whether we have to resolve any
				   // indirect rows generated from resolving
				   // constrained dofs.
  unsigned int added_rows = 0;

				   // first add the indices in an unsorted
				   // way and only keep track of the
				   // constraints that appear. They are
				   // resolved in a second step.
  for (unsigned int i = 0; i<n_local_dofs; ++i)
    {
      if (is_constrained(local_dof_indices[i]) == false)
	{
	  global_rows.global_row(added_rows)  = local_dof_indices[i];
	  global_rows.local_row(added_rows++) = i;
	  continue;
	}
      global_rows.insert_constraint(i);
    }
  global_rows.sort();

  const unsigned int n_constrained_rows = n_local_dofs-added_rows;
  for (unsigned int i=0; i<n_constrained_rows; ++i)
    {
      const unsigned int local_row = global_rows.last_constrained_local_row();
      Assert (local_row < n_local_dofs, ExcIndexRange(local_row, 0, n_local_dofs));
      const unsigned int global_row = local_dof_indices[local_row];
      Assert (is_constrained(global_row), ExcInternalError());
      const ConstraintLine & position =
	lines[lines_cache[calculate_line_index(global_row)]];
      if (position.inhomogeneity != 0)
	global_rows.set_last_row_inhomogeneous();
      else
	global_rows.pop_back();
      for (unsigned int q=0; q<position.entries.size(); ++q)
	global_rows.insert_index (position.entries[q].first,
				  local_row,
				  position.entries[q].second);

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
      const typename MatrixType::value_type new_diagonal
	= (std::fabs(local_matrix(local_row,local_row)) != 0 ?
	   std::fabs(local_matrix(local_row,local_row)) : average_diagonal);
      global_matrix.add(global_row, global_row, new_diagonal);
    }
}



template <typename SparsityType>
inline
void
ConstraintMatrix::
  make_sorted_row_list (const std::vector<unsigned int> &local_dof_indices,
			const bool                       keep_constrained_entries,
			SparsityType                    &sparsity_pattern,
			std::vector<unsigned int>       &actual_dof_indices) const
{
  const unsigned int n_local_dofs = local_dof_indices.size();
  unsigned int added_rows = 0;
  for (unsigned int i = 0; i<n_local_dofs; ++i)
    {
      if (is_constrained(local_dof_indices[i]) == false)
	{
	  actual_dof_indices[added_rows++] = local_dof_indices[i];
	  continue;
	}

      actual_dof_indices[n_local_dofs-i+added_rows-1] = i;
    }
  std::sort (actual_dof_indices.begin(), actual_dof_indices.begin()+added_rows);

  const unsigned int n_constrained_dofs = n_local_dofs-added_rows;
  for (unsigned int i=n_constrained_dofs; i>0; --i)
    {
      const unsigned int local_row = actual_dof_indices.back();
      actual_dof_indices.pop_back();
      const unsigned int global_row = local_dof_indices[local_row];
      const ConstraintLine & position =
	lines[lines_cache[calculate_line_index(global_row)]];
      for (unsigned int q=0; q<position.entries.size(); ++q)
	{
	  const unsigned int new_index = position.entries[q].first;
	  if (actual_dof_indices[actual_dof_indices.size()-i] < new_index)
	    actual_dof_indices.insert(actual_dof_indices.end()-i+1,new_index);
	  else
	    {
	      std::vector<unsigned int>::iterator it =
		std::lower_bound(actual_dof_indices.begin(),
				 actual_dof_indices.end()-i+1,
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
}



template <typename SparsityType>
inline
void
ConstraintMatrix::
  make_sorted_row_list (const Table<2,bool>                &dof_mask,
			const std::vector<unsigned int>    &local_dof_indices,
			const bool                          keep_constrained_entries,
			SparsityType                       &sparsity_pattern,
			internals::GlobalRowsFromLocal     &global_rows) const
{
				   // cache whether we have to resolve any
				   // indirect rows generated from resolving
				   // constrained dofs.
  unsigned int added_rows = 0;
  const unsigned int n_local_dofs = local_dof_indices.size();

  for (unsigned int i = 0; i<n_local_dofs; ++i)
    {
      if (is_constrained(local_dof_indices[i]) == false)
	{
	  global_rows.global_row(added_rows)  = local_dof_indices[i];
	  global_rows.local_row(added_rows++) = i;
	  continue;
	}
      global_rows.insert_constraint(i);
    }
  global_rows.sort();

  const unsigned int n_constrained_rows = n_local_dofs-added_rows;
  for (unsigned int i=0; i<n_constrained_rows; ++i)
    {
      const unsigned int local_row = global_rows.last_constrained_local_row();
      Assert (local_row < n_local_dofs, ExcIndexRange(local_row, 0, n_local_dofs));
      const unsigned int global_row = local_dof_indices[local_row];
      global_rows.pop_back();
      const ConstraintLine & position =
	lines[lines_cache[calculate_line_index(global_row)]];
      for (unsigned int q=0; q<position.entries.size(); ++q)
	global_rows.insert_index (position.entries[q].first,
				  local_row,
				  position.entries[q].second);

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
}



				// Resolve the constraints from the vector and
				// apply inhomogeneities.
inline
double
ConstraintMatrix::
  resolve_vector_entry (const unsigned int                    i,
			const internals::GlobalRowsFromLocal &global_rows,
			const Vector<double>                 &local_vector,
			const std::vector<unsigned int>      &local_dof_indices,
			const FullMatrix<double>             &local_matrix) const
{
  const unsigned int loc_row = global_rows.local_row(i);
  const unsigned int n_inhomogeneous_rows = global_rows.n_inhomogeneous_rows;
  double val = 0;
				// has a direct contribution from some local
				// entry. If we have inhomogeneous
				// constraints, compute the contribution of
				// the inhomogeneity in the current row.
  if (loc_row != numbers::invalid_unsigned_int)
    {
      val = local_vector(loc_row);
      for (unsigned int i=0; i<n_inhomogeneous_rows; ++i)
	val -= (lines[lines_cache[calculate_line_index(local_dof_indices
						       [global_rows.inhomogeneity(i)])]].
		inhomogeneity *
		local_matrix(loc_row, global_rows.inhomogeneity(i)));
    }

				// go through the indirect contributions
  for (unsigned int q=0; q<global_rows.size(i); ++q)
    {
      const unsigned int loc_row_q = global_rows.local_row(i,q);
      double add_this = local_vector (loc_row_q);
      for (unsigned int k=0; k<n_inhomogeneous_rows; ++k)
	add_this -= (lines[lines_cache[calculate_line_index(local_dof_indices
							    [global_rows.inhomogeneity(k)])]].
		     inhomogeneity *
		     local_matrix(loc_row_q,global_rows.inhomogeneity(k)));
      val += add_this * global_rows.constraint_value(i,q);
    }
  return val;
}


				   // internal implementation for
				   // distribute_local_to_global for
				   // standard (non-block) matrices
template <typename MatrixType, typename VectorType>
void
ConstraintMatrix::distribute_local_to_global (
  const FullMatrix<double>        &local_matrix,
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
  typedef typename MatrixType::value_type number;
  const bool use_dealii_matrix =
    types_are_equal<MatrixType,SparseMatrix<number> >::value;

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
  internals::GlobalRowsFromLocal global_rows (n_local_dofs);
  make_sorted_row_list (local_matrix, local_dof_indices, global_matrix,
			global_rows);

  const unsigned int n_actual_dofs = global_rows.size();

				   // create arrays for the column data
				   // (indices and values) that will then be
				   // written into the matrix. Shortcut for
				   // deal.II sparse matrix
  std::vector<unsigned int> cols;
  std::vector<number>       vals;
  SparseMatrix<number> * sparse_matrix
    = dynamic_cast<SparseMatrix<number> *>(&global_matrix);
  if (use_dealii_matrix == false)
    {
      cols.resize (n_actual_dofs);
      vals.resize (n_actual_dofs);
    }
  else
    Assert (sparse_matrix != 0, ExcInternalError());

				   // now do the actual job. go through all
				   // the global rows that we will touch and
				   // call resolve_matrix_row for each of
				   // those.
  for (unsigned int i=0; i<n_actual_dofs; ++i)
    {
      const unsigned int row = global_rows.global_row(i);

				   // calculate all the data that will be
				   // written into the matrix row.
      if (use_dealii_matrix == false)
	{
	  unsigned int * col_ptr = &cols[0];
	  number * val_ptr = &vals[0];
	  resolve_matrix_row (global_rows, i, 0, n_actual_dofs,
			      local_matrix, col_ptr, val_ptr);
	  const unsigned int n_values = col_ptr - &cols[0];
	  Assert (n_values == (unsigned int)(val_ptr - &vals[0]),
		  ExcInternalError());
	  if (n_values > 0)
	    global_matrix.add(row, n_values, &cols[0], &vals[0], false, true);
	}
      else
	resolve_matrix_row (global_rows, i, 0, n_actual_dofs,
			    local_matrix, sparse_matrix);

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
	  const double val = resolve_vector_entry (i, global_rows,
						   local_vector,
						   local_dof_indices,
						   local_matrix);

	  if (val != 0)
	    global_vector(row) += static_cast<typename VectorType::value_type>(val);
	}
    }
}



				   // similar function as above, but now
				   // specialized for block matrices. See
				   // the other function for additional
				   // comments.
template <typename MatrixType, typename VectorType>
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<double>        &local_matrix,
			    const Vector<double>            &local_vector,
                            const std::vector<unsigned int> &local_dof_indices,
                            MatrixType                      &global_matrix,
			    VectorType                      &global_vector,
			    internal::bool2type<true>) const
{
  const bool use_vectors = (local_vector.size() == 0 &&
			    global_vector.size() == 0) ? false : true;
  typedef typename MatrixType::value_type number;
  const bool use_dealii_matrix =
    types_are_equal<MatrixType,BlockSparseMatrix<number> >::value;

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
  internals::GlobalRowsFromLocal global_rows (n_local_dofs);
  make_sorted_row_list (local_matrix, local_dof_indices, global_matrix,
			global_rows);
  const unsigned int n_actual_dofs = global_rows.size();

  std::vector<unsigned int> global_indices;
  if (use_vectors == true)
    {
      global_indices.resize(n_actual_dofs);
      for (unsigned int i=0; i<n_actual_dofs; ++i)
	global_indices[i] = global_rows.global_row(i);
    }

				   // additional construct that also takes
				   // care of block indices.
  const unsigned int num_blocks   = global_matrix.n_block_rows();
  std::vector<unsigned int> block_starts(num_blocks+1, n_actual_dofs);
  internals::make_block_starts (global_matrix, global_rows, block_starts);

  std::vector<unsigned int> cols;
  std::vector<number>       vals;
  if (use_dealii_matrix == false)
    {
      cols.resize (n_actual_dofs);
      vals.resize (n_actual_dofs);
    }

				   // the basic difference to the non-block
				   // variant from now onwards is that we go
				   // through the blocks of the matrix
				   // separately, which allows us to set the
				   // block entries individually
  for (unsigned int block=0; block<num_blocks; ++block)
    {
      const unsigned int next_block = block_starts[block+1];
      for (unsigned int i=block_starts[block]; i<next_block; ++i)
	{
	  const unsigned int row = global_rows.global_row(i);

	  for (unsigned int block_col=0; block_col<num_blocks; ++block_col)
	    {
	      const unsigned int start_block = block_starts[block_col],
		end_block = block_starts[block_col+1];
	      if (use_dealii_matrix == false)
		{
		  unsigned int * col_ptr = &cols[0];
		  number * val_ptr = &vals[0];
		  resolve_matrix_row (global_rows, i, start_block,
				      end_block, local_matrix, col_ptr, val_ptr);
		  const unsigned int n_values = col_ptr - &cols[0];
		  Assert (n_values == (unsigned int)(val_ptr - &vals[0]),
			  ExcInternalError());
		  if (n_values > 0)
		    global_matrix.block(block, block_col).add(row, n_values,
							      &cols[0], &vals[0],
							      false, true);
		}
	      else
		{
		  SparseMatrix<number> * sparse_matrix
		    = dynamic_cast<SparseMatrix<number> *>(&global_matrix.block(block,
										block_col));
		  Assert (sparse_matrix != 0, ExcInternalError());
		  resolve_matrix_row (global_rows, i, start_block,
				      end_block, local_matrix, sparse_matrix);
		}
	    }

	  if (use_vectors == true)
	    {
	      const double val = resolve_vector_entry (i, global_rows,
						       local_vector,
						       local_dof_indices,
						       local_matrix);

	      if (val != 0)
		global_vector(global_indices[i]) +=
		  static_cast<typename VectorType::value_type>(val);
	    }
	}
    }
}



template <typename SparsityType>
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
      make_sorted_row_list (local_dof_indices, keep_constrained_entries,
			    sparsity_pattern, actual_dof_indices);
      const unsigned int n_actual_dofs = actual_dof_indices.size();

				   // now add the indices we collected above
				   // to the sparsity pattern. Very easy
				   // here - just add the same array to all
				   // the rows...
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
  internals::GlobalRowsFromLocal global_rows (n_local_dofs);
  make_sorted_row_list (dof_mask, local_dof_indices, keep_constrained_entries,
			sparsity_pattern, global_rows);
  const unsigned int n_actual_dofs = global_rows.size();

				   // create arrays for the column indices
				   // that will then be written into the
				   // sparsity pattern.
  std::vector<unsigned int> cols (n_actual_dofs);

  for (unsigned int i=0; i<n_actual_dofs; ++i)
    {
      std::vector<unsigned int>::iterator col_ptr = cols.begin();
      const unsigned int row = global_rows.global_row(i);
      resolve_matrix_row (global_rows, i, 0, n_actual_dofs,
			  dof_mask, col_ptr);

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

  if (dof_mask_is_active == false)
    {
      std::vector<unsigned int> actual_dof_indices (n_local_dofs);
      make_sorted_row_list (local_dof_indices, keep_constrained_entries,
			    sparsity_pattern, actual_dof_indices);
      const unsigned int n_actual_dofs = actual_dof_indices.size();

				   // additional construct that also takes
				   // care of block indices.
      std::vector<unsigned int> block_starts(num_blocks+1, n_actual_dofs);
      internals::make_block_starts (sparsity_pattern, actual_dof_indices,
				    block_starts);

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
  internals::GlobalRowsFromLocal global_rows (n_local_dofs);
  make_sorted_row_list (dof_mask, local_dof_indices, keep_constrained_entries,
			sparsity_pattern, global_rows);
  const unsigned int n_actual_dofs = global_rows.size();

				   // additional construct that also takes
				   // care of block indices.
  std::vector<unsigned int> block_starts(num_blocks+1, n_actual_dofs);
  internals::make_block_starts(sparsity_pattern, global_rows,
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
	  const unsigned int row = global_rows.global_row(i);
	  for (unsigned int block_col=0; block_col<num_blocks; ++block_col)
	    {
	      const unsigned int begin_block = block_starts[block_col],
		end_block = block_starts[block_col+1];
	      std::vector<unsigned int>::iterator col_ptr = cols.begin();
	      resolve_matrix_row (global_rows, i, begin_block, end_block,
				  dof_mask, col_ptr);

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
