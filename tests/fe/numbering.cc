// $Id$
// Author: Wolfgang Bangerth, 2001
//
// Check the numbering of continuous Lagrange finite elements. it
// constructs and independent numbering and compares it with the
// result of two functions from the library

#include <base/logstream.h>
#include <fe/fe_q.h>
#include <fe/fe_tools.h>
#include <vector>
#include <fstream>


std::ofstream logfile ("numbering.output");


template <int dim>
void check (const FE_Q<dim> &fe)
{
  Assert (fe.n_components() == 1, ExcInternalError());
  
  std::vector<unsigned int> hierarchic_to_lexicographic_numbering (fe.dofs_per_cell);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
				   // polynomial degree
  const unsigned int degree = fe.dofs_per_line+1;
				   // number of grid points in each direction
  const unsigned int n = degree+1;

  switch (dim)
    {
      case 1:
      {
	hierarchic_to_lexicographic_numbering[0] = 0;
	hierarchic_to_lexicographic_numbering[1] = dofs_per_cell-1;
	for (unsigned int i=2; i<dofs_per_cell; ++i)
	  hierarchic_to_lexicographic_numbering[i] = i-1;

	break;
      };

      case 2:
      {
	unsigned int next_index = 0;
					 // first the four vertices
	hierarchic_to_lexicographic_numbering[next_index++] = 0;
	hierarchic_to_lexicographic_numbering[next_index++] = n-1;
	hierarchic_to_lexicographic_numbering[next_index++] = n*n-1;
	hierarchic_to_lexicographic_numbering[next_index++] = n*(n-1);
					 // first line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = 1+i;
					 // second line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = (2+i)*n-1;
					 // third line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = n*(n-1)+i+1;
					 // fourth line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = (1+i)*n;
					 // inside quad
	Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
		ExcInternalError());
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    hierarchic_to_lexicographic_numbering[next_index++] = n*(i+1)+j+1;

	break;
      };

      case 3:
      {
	unsigned int next_index = 0;
					 // first the eight vertices
	hierarchic_to_lexicographic_numbering[next_index++] = 0;
	hierarchic_to_lexicographic_numbering[next_index++] = n-1;
	hierarchic_to_lexicographic_numbering[next_index++] = (n-1)*(n*n+1);
	hierarchic_to_lexicographic_numbering[next_index++] = (n-1)*n*n;
	hierarchic_to_lexicographic_numbering[next_index++] = n*(n-1);
	hierarchic_to_lexicographic_numbering[next_index++] = n*n-1;
	hierarchic_to_lexicographic_numbering[next_index++] = n*n*n-1;
	hierarchic_to_lexicographic_numbering[next_index++] = (n-1)*(n*n+n);

	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = 1+i;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = n-1+(i+1)*n*n;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = n*n*(n-1)+i+1;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n;

	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = 1+i+n*(n-1);
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = n-1+(i+1)*n*n+n*(n-1);
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = n*n*(n-1)+i+1+n*(n-1);
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n+n*(n-1);

	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = n-1+(i+1)*n;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = (n-1)*(n*n+1)+(i+1)*n;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  hierarchic_to_lexicographic_numbering[next_index++] = (n-1)*n*n+(i+1)*n;

					 // inside quads
	Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
		ExcInternalError());
					 // quad 1
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n+j+1;
					 // quad 2
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n+n*(n-1)+j+1;
					 // quad 3
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    hierarchic_to_lexicographic_numbering[next_index++] = n*(i+1)+j+1;
					 // quad 4
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n+n-1+n*(j+1);
					 // quad 5
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    hierarchic_to_lexicographic_numbering[next_index++] = (n-1)*n*n+n*(i+1)+j+1;
					 // quad 6
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n+n*(j+1);

					 // inside hex
	Assert (fe.dofs_per_hex == fe.dofs_per_quad*fe.dofs_per_line,
		ExcInternalError());
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    for (unsigned int k=0; k<fe.dofs_per_line; ++k)
	      hierarchic_to_lexicographic_numbering[next_index++]
		= n*n*(i+1)+n*(j+1)+k+1;
	
	break;
      };
       

      default:
	    Assert (false, ExcNotImplemented());
    };
  
				   // now check with data from the
				   // lib: there we have the mapping
				   // the other way round, and we
				   // check that the concatenation of
				   // the two mappings is the
				   // identity. output the two maps to
				   // generate some output for
				   // automatic comparison
  std::vector<unsigned int> l2h (fe.dofs_per_cell);
  FETools::lexicographic_to_hierarchic_numbering (fe, l2h);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    {
      Assert (l2h[hierarchic_to_lexicographic_numbering[i]] == i,
	      ExcInternalError());
      logfile << dim << "d, degree=" << degree << ": "
	      << l2h[i]
	      << ' '
	      << hierarchic_to_lexicographic_numbering[i]
	      << std::endl;
    };

				   // finally, we also have the
				   // forward map in the lib, so check
				   // for equality
  std::vector<unsigned int> h2l (fe.dofs_per_cell);
  FETools::hierarchic_to_lexicographic_numbering (fe, h2l);
  Assert (hierarchic_to_lexicographic_numbering == h2l,
	  ExcInternalError());
};



template <int dim>
void check_dim ()
{
  for (unsigned int degree=1; degree<6; ++degree)
    {
      FE_Q<dim> fe(degree);
      check (fe);
    };
};


int main ()
{
  logfile.precision (2);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console (0);

  check_dim<1> ();
  check_dim<2> ();
  check_dim<3> ();
};


