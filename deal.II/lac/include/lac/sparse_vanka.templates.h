// $Id$
// Copyright Guido Kanschat, 1999

#include <lac/sparse_vanka.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>

#include <map>



template<typename number>
SparseVanka<number>::SparseVanka(const SparseMatrix<number> &M,
				 const vector<bool>         &selected)
		:
		matrix(&M),
		selected(selected),
		conserve_mem(false),
		inverses(M.m(),0)
{
  Assert (M.m() == M.n(),
	  ExcMatrixNotSquare ());
}


template<typename number>
SparseVanka<number>::~SparseVanka()
{
  vector<SmartPointer<FullMatrix<float> > >::iterator i;
  for(i=inverses.begin();i!=inverses.end();++i)
    {
      FullMatrix<float>* p = *i;
      *i = 0;
      if (p != 0) delete p;
    }
}



template<typename number>
template<typename number2>
void
SparseVanka<number>::operator ()(Vector<number2>       &dst,
				 const Vector<number2> &src) const
{
				   // first set output vector to zero
  dst.clear ();
				   // first define an alias to the sparsity
				   // pattern of the matrix, since this
				   // will be used quite often
  const SparseMatrixStruct &structure
    = matrix->get_sparsity_pattern();
  
				   // space to be used for local
				   // systems. allocate as much memory
				   // as is the maximum. this
				   // eliminates the need to
				   // re-allocate memory inside the
				   // loop.
  FullMatrix<float> local_matrix (structure.max_row_length(),
				  structure.max_row_length());
  Vector<float> b (structure.max_row_length());
  Vector<float> x (structure.max_row_length());
  
				   // traverse all rows of the matrix
  for (unsigned int row=0; row< matrix->m() ; ++row)
    {
				       // but skip those that were not selected
      if (!selected[row])
	continue;

				       // shall we reconstruct the system
				       // of equations for this DoF?
      bool build_matrix = true;
      const unsigned int row_length = structure.row_length(row);
      
      if (!conserve_mem)
	{
	  if (inverses[row] == 0)
					     // inverse not yet built
	    {
	      FullMatrix<float> *p = new FullMatrix<float>;
	      inverses[row] = p;
	    } else {
					       // inverse already built
	      build_matrix = false;
	    }
	}

				       // if we don't store the
				       // inverse matrices, then alias
				       // the entry in the global
				       // vector to the local matrix
				       // to be used
      if (conserve_mem == true)
	{
	  inverses[row] = &local_matrix;
	  inverses[row]->reinit (row_length, row_length);
	}
      else
	if (build_matrix)
	  inverses[row]->reinit (row_length, row_length);
      
      b.reinit (row_length);
      x.reinit (row_length);

				       // mapping between:
				       // 1 column number of all
				       //   entries in this row, and
				       // 2 the position within this
				       //   row (as stored in the
				       //   sparsematrixstruct object
				       //
				       // since we do not explicitely
				       // consider nonsysmmetric sparsity
				       // patterns, the first element
				       // of each entry simply denotes
				       // all degrees of freedom that
				       // couple with #row#.
      map<unsigned int, unsigned int> local_index;
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
					   // row corresponding to DoF #row#.
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
					       // (that couples with #irow#,
					       // which itself couples with
					       // #row#) also couples with
					       // #row#.
	      const map<unsigned int, unsigned int>::const_iterator js
		= local_index.find(col);
					       // if not, then still use
					       // this dof to modify the rhs
					       //
					       // note that if so, we already
					       // have copied the entry above
	      if (js == local_index.end())
		{
		  b(i) -= matrix->raw_entry(irow,j) * dst(col);
		} else {
						   // if so, then build the
						   // matrix out of it
		  if (build_matrix)
		    (*inverses[row])(i,js->second) = matrix->raw_entry(irow,j);
		}
	    }
	}
				       // Compute new values
      if (build_matrix)
	inverses[row]->gauss_jordan();
      inverses[row]->vmult(x,b);
      
				       // Distribute new values
      for (map<unsigned int, unsigned int>::const_iterator is=local_index.begin();
	   is!=local_index.end(); ++is)
	{
	  const unsigned int irow = is->first;
	  const unsigned int i = is->second;
	  dst(irow) = x(i);
	};

				       // if we don't store the
				       // inverses, then unalias the
				       // local matrix
      if (conserve_mem == true)
	inverses[row] = 0;
    }
}

