// $Id$
// Copyright Guido Kanschat, 1999

#include <lac/sparsevanka.h>
#include <lac/fullmatrix.h>

#include <map>

template<typename number>
SparseVanka<number>::SparseVanka(const SparseMatrix<number>& M,
				 const vector<int>& indices)
		:
		matrix(&M), indices(indices)
{}

template<typename number>
template<typename number2>
void
SparseVanka<number>::apply(Vector<number2>& dst) const
{
  for (unsigned int global_i=0; global_i<indices.size() ; ++global_i)
    {
      unsigned int row = indices[global_i];
      const SparseMatrixStruct& structure = matrix->get_sparsity_pattern();
      unsigned int n = structure.row_length(row);
      
      FullMatrix<number> A(n);
      Vector<number> b(n);
      Vector<number> x(n);
      
      map<unsigned int, unsigned int> local_index;

				       // Build local index

      for (unsigned int i=0;i<n;++i)
	local_index.insert(pair<unsigned int, unsigned int>
			   (structure.column_number(row, i), i));

				       // Build local matrix

      for (map<unsigned int, unsigned int>::iterator is=local_index.begin();
	   is!=local_index.end();++is)
	{
	  unsigned int irow = is->first;
	  unsigned int i = is->second;
	  unsigned int n = structure.row_length(irow);
	  
	  b(i) = dst(irow);
	  
	  for (unsigned int j=0;j<n;++j)
	    {
	      unsigned int col = structure.column_number(irow, j);
	      map<unsigned int, unsigned int>::iterator js
		= local_index.find(col);
	      if (js == local_index.end())
		{
		  b(i) -= matrix->raw_entry(irow,col) * dst(col);
		} else {
		  A(i,j) = matrix->raw_entry(irow,col);
		}
	    }
	}
				       // Compute new values
      A.gauss_jordan();
      A.vmult(x,b);
      
				       // Distribute new values
      for (map<unsigned int, unsigned int>::iterator is=local_index.begin();
	   is!=local_index.end();++is)
	{
	  unsigned int irow = is->first;
	  unsigned int i = is->second;
	  dst(irow) = x(i);
	}     
    }
}

