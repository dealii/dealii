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
plot_transformation(FiniteElement<dim>& finel,
		    DoFHandler<dim>::cell_iterator& c,
		    const char* name)
{
  const unsigned int div = 20;

  QIteratedTrapez<dim> q(div);
  FEValues<dim> fe(finel, q,
		   UpdateFlags(update_q_points | update_JxW_values));

  fe.reinit(c);
  
  sprintf(fname, "%s.dat", name);
  ofstream gnuplot(fname);

  unsigned int k=0;
  for (unsigned int m=0;m<=div;++m)
    {
      for (unsigned int n=0;n<=div;++n)
	{
	  gnuplot << fe.quadrature_point(k);
	  double J = fe.JxW(k) / q.weight(k);
	  gnuplot << ' ' << J << endl;
	  k++;
	}
      gnuplot << endl;
    }     
}


main()
{
  Triangulation<2> tr;
  FEQ1<2> q1;
  DoFHandler<2> dof(&tr);

  GridGenerator::hyper_cube(tr, -1., 1.);
  dof.distribute_dofs(q1);
  
  DoFHandler<2>::cell_iterator c = dof.begin();

  Point<2>& v = c->vertex(2);
  
  v(0) = 3.;
  v(1) = 2.;
  
  plot_transformation(q1,c,"Transform-Q1");
}
