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
plot_shape_functions(Mapping<dim>& mapping,
		     FiniteElement<dim>& finel,
		     const char* name)
{
  Triangulation<dim> tr;
  DoFHandler<dim> dof(tr);
  GridGenerator::hyper_cube(tr, 0., 1.);
  DoFHandler<dim>::cell_iterator c = dof.begin();
  dof.distribute_dofs(finel);

  const unsigned int div = 11;

  QTrapez<1> q_trapez;
  QIterated<dim> q(q_trapez, div);
  FEValues<dim> fe(mapping, finel, q, UpdateFlags(update_values));
  fe.reinit(c);
  
  sprintf(fname, "Shapes%dd-%s.output", dim, name);
  ofstream gnuplot(fname);
  gnuplot.setf(ios::fixed);
  gnuplot.precision (2);

  unsigned int k=0;
  for (unsigned int mz=0;mz<=((dim>2) ? div : 0) ;++mz)
    {
      for (unsigned int my=0;my<=((dim>1) ? div : 0) ;++my)
	{
	  for (unsigned int mx=0;mx<=div;++mx)
	    {
	      gnuplot << q.point(k);
	  
	      for (unsigned int i=0;i<finel.dofs_per_cell;++i)
		{
		  gnuplot << " " << fe.shape_value(i,k) + 1.;
		}
	      gnuplot << std::endl;
	      k++;
	    }
	  gnuplot << std::endl;
	}
      gnuplot << std::endl;
    }
}



template<int dim>
inline void
plot_face_shape_functions(Mapping<dim>& mapping,
			  FiniteElement<dim>& finel,
			  const char* name)
{
  Triangulation<dim> tr;
  DoFHandler<dim> dof(tr);
  GridGenerator::hyper_cube(tr, 0., 1.);
  tr.refine_global(1);
  DoFHandler<dim>::active_cell_iterator c = dof.begin_active();
  ++c;
  c->set_refine_flag();
  tr.execute_coarsening_and_refinement ();
  c = dof.begin_active();

  dof.distribute_dofs(finel);

  const unsigned int div = 4;

  QTrapez<1> q_trapez;
  QIterated<dim-1> q(q_trapez, div);
  FEFaceValues<dim> fe(mapping, finel, q, UpdateFlags(update_values
    | update_q_points));
  FESubfaceValues<dim> sub(mapping, finel, q, UpdateFlags(update_values
    | update_q_points));

  sprintf(fname, "ShapesFace%dd-%s.output", dim, name);
  ofstream gnuplot(fname);
  gnuplot.setf(ios::fixed);
  gnuplot.precision (2);
  
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    {
      if (!c->face(f)->has_children())
	{
	  fe.reinit(c, f);
  
	  unsigned int k=0;
	  for (unsigned int my=0;my<=((dim>2) ? div : 0) ;++my)
	    {
	      for (unsigned int mx=0;mx<=div;++mx)
		{
		  gnuplot << fe.quadrature_point(k);
		  
		  for (unsigned int i=0;i<finel.dofs_per_cell;++i)
		    {
		      gnuplot << " " << fe.shape_value(i,k) + 1.;
		    }
		  gnuplot << std::endl;
		  k++;
		}
	      gnuplot << std::endl;
	    }
	  gnuplot << std::endl;
	} else {
	  for (unsigned int s=0;s<GeometryInfo<dim>::subfaces_per_face; ++s)
	    {
	      sub.reinit(c, f, s);
	      
	      unsigned int k=0;
	      for (unsigned int my=0;my<=((dim>2) ? div : 0) ;++my)
		{
		  for (unsigned int mx=0;mx<=div;++mx)
		    {
		      gnuplot << sub.quadrature_point(k);
		      
		      for (unsigned int i=0;i<finel.dofs_per_cell;++i)
			{
			  gnuplot << " " << sub.shape_value(i,k) + 1.;
			}
		      gnuplot << std::endl;
		      k++;
		    }
		  gnuplot << std::endl;
		}
	      gnuplot << std::endl;	      
	    }
	}
    }
}


template<>
void plot_face_shape_functions (Mapping<1>&,
				FiniteElement<1>&,
				const char*)
{}


template<int dim>
void plot_FE_Q_shape_functions()
{
  MappingQ1<dim> m;
  FE_Q<dim> q1(1);
  plot_shape_functions(m, q1, "Q1");
  plot_face_shape_functions(m, q1, "Q1");
  FE_Q<dim> q2(2);
  plot_shape_functions(m, q2, "Q2");
  plot_face_shape_functions(m, q2, "Q2");
  FE_Q<dim> q3(3);
  plot_shape_functions(m, q3, "Q3");
  plot_face_shape_functions(m, q3, "Q3");
  FE_Q<dim> q4(4);
  plot_shape_functions(m, q4, "Q4");
  plot_face_shape_functions(m, q4, "Q4");
//    FE_Q<dim> q5(5);
//    plot_shape_functions(m, q5, "Q5");
//    FE_Q<dim> q6(6);
//    plot_shape_functions(m, q6, "Q6");
//    FE_Q<dim> q7(7);
//    plot_shape_functions(m, q7, "Q7");
//    FE_Q<dim> q8(8);
//    plot_shape_functions(m, q8, "Q8");
//    FE_Q<dim> q9(9);
//    plot_shape_functions(m, q9, "Q9");
//    FE_Q<dim> q10(10);
//    plot_shape_functions(m, q10, "Q10");
}


template<int dim>
void plot_FE_DGQ_shape_functions()
{
  MappingQ1<dim> m;
  FE_DGQ<dim> q1(1);
  plot_shape_functions(m, q1, "DGQ1");
  plot_face_shape_functions(m, q1, "DGQ1");
  FE_DGQ<dim> q2(2);
  plot_shape_functions(m, q2, "DGQ2");
  plot_face_shape_functions(m, q2, "DGQ2");
  FE_DGQ<dim> q3(3);
  plot_shape_functions(m, q3, "DGQ3");
  plot_face_shape_functions(m, q3, "DGQ3");
  FE_DGQ<dim> q4(4);
  plot_shape_functions(m, q4, "DGQ4");
  plot_face_shape_functions(m, q4, "DGQ4");
//    FE_DGQ<dim> q5(5);
//    plot_shape_functions(m, q5, "DGQ5");
//    FE_DGQ<dim> q6(6);
//    plot_shape_functions(m, q6, "DGQ6");
//    FE_DGQ<dim> q7(7);
//    plot_shape_functions(m, q7, "DGQ7");
//    FE_DGQ<dim> q8(8);
//    plot_shape_functions(m, q8, "DGQ8");
//    FE_DGQ<dim> q9(9);
//    plot_shape_functions(m, q9, "DGQ9");
//    FE_DGQ<dim> q10(10);
//    plot_shape_functions(m, q10, "DGQ10");
}


int
main()
{
  ofstream logfile ("shapes.output");
  logfile.precision (2);
  logfile.setf(ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  plot_FE_Q_shape_functions<1>();
  plot_FE_Q_shape_functions<2>();
  plot_FE_DGQ_shape_functions<2>();
//  plot_FE_Q_shape_functions<3>();

				   // FESystem test.
  MappingQ1<2> m;
  FESystem<2> q2_q3(FE_Q<2>(2), 1,
		    FE_Q<2>(3), 1);
//  plot_shape_functions(m, q2_q3, "Q2_Q3");
  
  return 0;
}



