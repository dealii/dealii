// $Id$
// (c) Guido Kanschat
//

#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/full_matrix.templates.h>
#include <lac/solver_richardson.h>
#include <lac/vector_memory.h>
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
#include <iomanip>
#include <fstream>
#include <strstream>
#include <string>


template <int dim>
void compute_embedding (Triangulation<dim>& tr_coarse,
			Triangulation<dim>& tr_fine,
			const FiniteElement<dim>& fe_coarse,
			const FiniteElement<dim>& fe_fine,
			const char* name)
{
  cerr << name << '<' << dim << '>' << endl;
  
  DoFHandler<dim> dof_coarse(tr_coarse);
  dof_coarse.distribute_dofs(fe_coarse);
  DoFHandler<dim> dof_fine(tr_fine);
  dof_fine.distribute_dofs(fe_fine);
  
  FullMatrix<double> A(dof_fine.n_dofs());
  Vector<long double> f(dof_fine.n_dofs());
  Vector<long double> u(dof_fine.n_dofs());
  vector<vector<vector<double> > >
    result(2<<dim,
	   vector<vector<double> >(fe_coarse.dofs_per_cell,
				   vector<double>(fe_fine.dofs_per_cell, 0.)));

  DoFHandler<dim>::active_cell_iterator coarse_cell
    = dof_coarse.begin_active();
  vector<unsigned int> indices(fe_fine.dofs_per_cell);

  MappingQ1<dim> mapping;
  QGauss<dim> q_fine(12);
  
  FEValues<dim> fine (mapping, fe_fine, q_fine,
		      update_q_points
		      | update_JxW_values
		      | update_values);


  DoFHandler<dim>::active_cell_iterator c;

  for (unsigned int coarse_no=0;coarse_no<fe_coarse.dofs_per_cell;
       ++coarse_no)
    {
      A.clear();
      f.clear();
      u.clear();

      for (c = dof_fine.begin_active();
	   c != dof_fine.end();
	   ++c)
	{
	  c->get_dof_indices(indices);
	  fine.reinit(c);
					   // Build mass matrix and RHS
	  
	  Quadrature<dim> q_coarse (fine.get_quadrature_points(),
				    fine.get_JxW_values());
	  FEValues<dim> coarse (mapping, fe_coarse, q_coarse,
				update_values);
	  coarse.reinit(coarse_cell);
	  
	  for (unsigned int k=0;k<fine.n_quadrature_points;++k)
	    {
	      double dx = fine.JxW(k);
	      for (unsigned int i=0;i<fe_fine.dofs_per_cell;++i)
		{
		  double v = fine.shape_value(i,k);
		  f(indices[i]) += dx *
		    v * coarse.shape_value(coarse_no, k);
		  for (unsigned int j=0;j<fe_fine.dofs_per_cell;++j)
		    A(indices[i],indices[j]) += dx*v*fine.shape_value(j,k);
		}
	    }
	}
//      A.print_formatted(cout, 2, false, 4, "~", 9);
      FullMatrix<double> P(A.n());
      P.invert(A);
      SolverControl control (100, 1.e-20, true, false);
      PrimitiveVectorMemory<Vector<long double> > mem;
      SolverRichardson<Vector<long double> > solver(control, mem);
      solver.solve(A,u,f,P);
      
      unsigned int cell_i=0;
      for (c = dof_fine.begin_active();
	   c != dof_fine.end();
	   ++c, ++cell_i)
	{
	  c->get_dof_indices(indices);
	  for (unsigned int i=0;i<fe_fine.dofs_per_cell;++i)
	    result[cell_i][coarse_no][i] = (fabs(u(indices[i])) > 1.e-16)
	      ? (27.*u(indices[i])) : 0.;
	}
    }

  for (unsigned int cell=0;cell<tr_fine.n_active_cells();++cell)
    {
      cout << "static double[] "
	   << name << "_"
	   << dim << "d_"
	   << cell << " =" << endl << '{' << endl;
      for (unsigned int i=0;i<fe_fine.dofs_per_cell;++i)
	{
	  for (unsigned int j=0;j<fe_coarse.dofs_per_cell;++j)
	    {
	      cout << ' ' << setprecision(10) << result[cell][j][i] << "/27.,";
	    }
	  cout << endl;
	}
      cout << "};" << endl << endl;
    }
}

#define PUSH_Q(k) FE_Q<dim> q ## k (k);\
 elements.push_back (ElPair(&q ## k, "q" # k))

#define PUSH_DGQ(k) FE_DGQ<dim> dgq ## k(k);\
 elements.push_back (ElPair(&dgq ## k, "dgq" # k))

template <int dim>
void loop ()
{
  Triangulation<dim> tr_coarse;
  Triangulation<dim> tr_fine;
  GridGenerator::hyper_cube (tr_coarse, 0, 1);
  GridGenerator::hyper_cube (tr_fine, 0, 1);
  tr_fine.refine_global(1);

  typedef pair<const FiniteElement<dim>*, const char*> ElPair ;
  vector <ElPair> elements(0);

        PUSH_DGQ(0);
        PUSH_DGQ(1);
        PUSH_DGQ(2);
//        PUSH_DGQ(3);
        PUSH_DGQ(5);
        PUSH_DGQ(6);
        PUSH_DGQ(7);
      //  PUSH_Q(0);
      PUSH_Q(1);
      PUSH_Q(2);
      PUSH_Q(3);
      PUSH_Q(4);

  char* name = new char[100];

				   // Embed all lower spaces into higher
  unsigned int n = elements.size();
  for (unsigned int i=0;i<n;++i)
    for (unsigned int j=i;j<=i;++j)
      {
	ostrstream os (name, 99);
	os << elements[i].second << "_into_"
	   << elements[j].second << "_refined" << ends;
	
	compute_embedding (tr_coarse, tr_fine,
			   *elements[i].first,
			   *elements[j].first,
			   name);
      }

  delete [] name;
}

int main ()
{
  loop<1> ();
  loop<2> ();
  loop<3> ();
}
