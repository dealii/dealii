// $Id$
// Copyright Guido Kanschat, 1999

#include <lac/sparse_vanka.h>
#include <lac/fullmatrix.h>

#include <map>

template<typename number>
SparseVanka<number>::SparseVanka(const SparseMatrix<number>& M,
				 const bit_vector& selected)
		:
		matrix(&M), selected(selected)
{}

template<typename number>
template<typename number2>
void
SparseVanka<number>::forward(Vector<number2>& dst,
			   const Vector<number2>& src) const
{
  for (unsigned int row=0; row< matrix->m() ; ++row)
    {
      if (!selected[row])
	continue;
      
      const SparseMatrixStruct& structure = matrix->get_sparsity_pattern();
      unsigned int n = structure.row_length(row);
      
      FullMatrix<float> A(n);
      Vector<float> b(n);
      Vector<float> x(n);
      
      map<unsigned int, unsigned int> local_index;

				       // Build local index

      for (unsigned int i=0;i<n;++i)
	local_index.insert(pair<unsigned int, unsigned int>
			   (structure.column_number(row, i), i));

//       for (map<unsigned int, unsigned int>::iterator is=local_index.begin();
// 	   is!=local_index.end();++is)
// 	cerr << "map " << is->first << '\t' << is->second << endl;
      
				       // Build local matrix

      for (map<unsigned int, unsigned int>::iterator is=local_index.begin();
	   is!=local_index.end();++is)
	{
	  unsigned int irow = is->first;
	  unsigned int i = is->second;
	  unsigned int n = structure.row_length(irow);
	  
	  b(i) = src(irow);
	  
	  for (unsigned int j=0;j<n;++j)
	    {
	      unsigned int col = structure.column_number(irow, j);
	      map<unsigned int, unsigned int>::iterator js
		= local_index.find(col);
	      if (js == local_index.end())
		{
		  b(i) -= matrix->raw_entry(irow,j) * dst(col);
		} else {
		  A(i,js->second) = matrix->raw_entry(irow,j);
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

