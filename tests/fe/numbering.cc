// $Id$
// Author: Wolfgang Bangerth, 2001
//
// Check the numbering of finite elements

#include <base/logstream.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
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
	Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line, ExcInternalError());
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    hierarchic_to_lexicographic_numbering[next_index++] = n*(i+1)+j+1;

	break;
      };

      default:
	    Assert (false, ExcNotImplemented());
    };
  
				   // now check with data from the lib
  const std::vector<unsigned int> &
    lexicographic_to_hierarchical_numbering = fe.get_renumber();
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    {
      Assert (lexicographic_to_hierarchical_numbering
	      [hierarchic_to_lexicographic_numbering[i]] == i,
	      ExcInternalError());
      logfile << dim << "d, degree=" << degree << ": "
	      << lexicographic_to_hierarchical_numbering[i]
	      << ' '
	      << hierarchic_to_lexicographic_numbering[i]
	      << std::endl;
    };
};



template <int dim>
void check_dim ()
{
  for (unsigned int degree=1; degree<4; ++degree)
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
};


