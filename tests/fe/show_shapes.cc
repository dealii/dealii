// $Id$
// (c) Guido Kanschat
//
// Show the shape functions implemented.

#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/fe_values.h>
#include <vector>
#include <fstream>
#include <string>

char fname[50];

template<int dim>
inline void
plot_shape_functions(FiniteElement<dim>& finel, const char* name)
{
  const unsigned int div = 20;

  QTrapez<1> q_trapez;
  QIterated<dim> q(q_trapez, div);
  FEValues<dim> fe(finel, q, UpdateFlags(update_values));

  sprintf(fname, "%s.dat", name);
  ofstream gnuplot(fname);
  gnuplot.setf(ios::fixed);
  gnuplot.precision (3);

  unsigned int k=0;
  for (unsigned int m=0;m<=div;++m)
    {
      for (unsigned int n=0;n<=div;++n)
	{
	  gnuplot << q.point(k);
	  
	  for (unsigned int i=0;i<finel.dofs_per_cell;++i)
	    {
	      gnuplot << " "<< fe.shape_value(i,k);
	    }
	  gnuplot << endl;
	  k++;
	}
      gnuplot << endl;
    }     
}


int
main()
{
  FEQ1<2> q1;
  plot_shape_functions(q1,"Q1");
  FEQ2<2> q2;
  plot_shape_functions(q2,"Q2");
  FEQ3<2> q3;
  plot_shape_functions(q3,"Q3");
  FEQ4<2> q4;
  plot_shape_functions(q4,"Q4");

  return 0;
}
