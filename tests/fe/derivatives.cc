// $Id$
// (c) Guido Kanschat
//
// Show the shape functions implemented.

#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/grid_generator.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1.h>
#include <fe/fe_values.h>
#include <vector>
#include <fstream>
#include <string>

char fname[50];

template<int dim>
inline void
plot_derivatives(Mapping<dim>& mapping,
		     FiniteElement<dim>& finel,
		     const char* name)
{
  deallog.push (name);
  
  Triangulation<dim> tr;
  DoFHandler<dim> dof(tr);
  GridGenerator::hyper_cube(tr, 2., 5.);
  DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(finel);

  const unsigned int div = 1;

  QTrapez<dim> q;
//  QIterated<dim> q(q_trapez, div);
  FEValues<dim> fe(mapping, finel, q, UpdateFlags(update_gradients
						  | update_second_derivatives));
  fe.reinit(c);

  unsigned int k=0;
  for (unsigned int mz=0;mz<=((dim>2) ? div : 0) ;++mz)
    {
      for (unsigned int my=0;my<=((dim>1) ? div : 0) ;++my)
	{
	  for (unsigned int mx=0;mx<=div;++mx)
	    {
	      deallog << q.point(k) << std::endl;
	  
	      for (unsigned int i=0;i<finel.dofs_per_cell;++i)
		{
		  deallog << "\tGrad " << fe.shape_grad(i,k);
		  deallog << "\t2nd " << fe.shape_2nd_derivative(i,k);
		  deallog << std::endl;
		}
	      k++;
	    }
	}
    }
  deallog.pop ();
}



template<int dim>
void plot_FE_Q_shape_functions()
{
  MappingQ1<dim> m;
  FE_Q<dim> q1(1);
  plot_derivatives(m, q1, "Q1");
//  plot_face_shape_functions(m, q1, "Q1");
  FE_Q<dim> q2(2);
  plot_derivatives(m, q2, "Q2");
//  plot_face_shape_functions(m, q2, "Q2");
  FE_Q<dim> q3(3);
  plot_derivatives(m, q3, "Q3");
//  plot_face_shape_functions(m, q3, "Q3");
  FE_Q<dim> q4(4);
  plot_derivatives(m, q4, "Q4");
//  plot_face_shape_functions(m, q4, "Q4");
//    FE_Q<dim> q5(5);
//    plot_derivatives(m, q5, "Q5");
//    FE_Q<dim> q6(6);
//    plot_derivatives(m, q6, "Q6");
//    FE_Q<dim> q7(7);
//    plot_derivatives(m, q7, "Q7");
//    FE_Q<dim> q8(8);
//    plot_derivatives(m, q8, "Q8");
//    FE_Q<dim> q9(9);
//    plot_derivatives(m, q9, "Q9");
//    FE_Q<dim> q10(10);
//    plot_derivatives(m, q10, "Q10");
}


template<int dim>
void plot_FE_DGQ_shape_functions()
{
  MappingQ1<dim> m;
  FE_DGQ<dim> q1(1);
  plot_derivatives(m, q1, "DGQ1");
//  plot_face_shape_functions(m, q1, "DGQ1");
  FE_DGQ<dim> q2(2);
  plot_derivatives(m, q2, "DGQ2");
//  plot_face_shape_functions(m, q2, "DGQ2");
  FE_DGQ<dim> q3(3);
  plot_derivatives(m, q3, "DGQ3");
//  plot_face_shape_functions(m, q3, "DGQ3");
  FE_DGQ<dim> q4(4);
  plot_derivatives(m, q4, "DGQ4");
//  plot_face_shape_functions(m, q4, "DGQ4");
//    FE_DGQ<dim> q5(5);
//    plot_derivatives(m, q5, "DGQ5");
//    FE_DGQ<dim> q6(6);
//    plot_derivatives(m, q6, "DGQ6");
//    FE_DGQ<dim> q7(7);
//    plot_derivatives(m, q7, "DGQ7");
//    FE_DGQ<dim> q8(8);
//    plot_derivatives(m, q8, "DGQ8");
//    FE_DGQ<dim> q9(9);
//    plot_derivatives(m, q9, "DGQ9");
//    FE_DGQ<dim> q10(10);
//    plot_derivatives(m, q10, "DGQ10");
}


int
main()
{
  std::ofstream logfile ("derivatives.output");
  logfile.precision (2);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  deallog.push ("1d");
  plot_FE_Q_shape_functions<1>();
  deallog.pop ();
  
  deallog.push ("2d");
  plot_FE_Q_shape_functions<2>();
  plot_FE_DGQ_shape_functions<2>();
  deallog.pop ();
  
  deallog.push ("3d");
//  plot_FE_Q_shape_functions<3>();
  deallog.pop ();
  

				   // FESystem test.
  MappingQ1<2> m;
  FESystem<2> q2_q3(FE_Q<2>(2), 1,
		    FE_Q<2>(3), 1);
//  plot_derivatives(m, q2_q3, "Q2_Q3");
  
  return 0;
}



