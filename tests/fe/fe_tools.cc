//--------------------------------------------------------------------------
//    fe_tools.cc,v 1.1 2003/11/28 15:03:26 guido Exp
//    Version: 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------------

// Test the cell matrices generated in FETools and the local renumbering vector.

#include "../tests.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#include <base/logstream.h>

#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>

#include <fe/fe_tools.h>


template <typename number>
void print_formatted (const FullMatrix<number> &A,
		      const unsigned int        precision,
		      const unsigned int        width)
{
  for (unsigned int i=0; i<A.m(); ++i)
    {
      for (unsigned int j=0; j<A.n(); ++j)
	{
	  if (A(i,j) != 0)
	    deallog << std::setw(width) << std::setprecision(precision)
		    << A(i,j);
	  else
	    deallog << std::setw(width) << std::setprecision(precision)
		    << "~";
	  deallog << ' ';
	};
      deallog << std::endl;
    };
}



template<int dim>
void test_embedding (const FiniteElement<dim>& fe)
{  
  const unsigned int n = fe.dofs_per_cell;
  const unsigned int nc= GeometryInfo<dim>::children_per_cell;
  
  FullMatrix<double> P[nc];
  for (unsigned int i=0;i<nc;++i)
    P[i].reinit(n,n);
  
  FETools::compute_embedding_matrices(fe, P);
  
  for (unsigned int i=0;i<nc;++i)
    {
      deallog << fe.get_name() << " embedding " << i << std::endl;
      print_formatted(P[i], 3, 6);
    }
}

  
template<int dim>
void test_projection (const FiniteElement<dim>& fe1,
		      const FiniteElement<dim>& fe2,
		      std::ostream& out)
{
  out << fe1.get_name() << " -> "
      << fe2.get_name() << std::endl;
  
  const unsigned int n1 = fe1.dofs_per_cell;
  const unsigned int n2 = fe2.dofs_per_cell;

  FullMatrix<double> P(n2,n1);

  FETools::get_projection_matrix(fe1, fe2, P);
  P.print_formatted(out, 2, false, 5);
}

  
template<int dim>
void test_projection (std::ostream& out)
{
  FE_DGQ<dim> q0(0);
  FE_DGQ<dim> q1(1);
  FE_DGQ<dim> q2(2);
  FE_DGQ<dim> q3(3);
  FE_DGQ<dim> q4(4);

  FE_DGP<dim> p0(0);
  FE_DGP<dim> p1(1);
  FE_DGP<dim> p2(2);
  FE_DGP<dim> p3(3);
  FE_DGP<dim> p4(4);

  test_embedding(q0);
  test_embedding(q1);
  test_embedding(q2);
  test_embedding(q3);
  test_embedding(p1);
  test_embedding(p2);
  test_embedding(p3);
  
  test_projection(p1,p0, out);
  test_projection(p0,p1, out);
  test_projection(p2,p1, out);
  test_projection(p1,p2, out);
  test_projection(p2,p0, out);
  test_projection(p0,p2, out);
  test_projection(p3,p2, out);
  test_projection(p2,p3, out);
  test_projection(p4,p3, out);
  test_projection(p3,p4, out);

  test_projection(q1,q0, out);
  test_projection(q2,q0, out);
  test_projection(q3,q0, out);
  test_projection(q4,q0, out);
  test_projection(q2,q1, out);
  test_projection(q1,q2, out);
  test_projection(q3,q2, out);
  test_projection(q4,q3, out);
}


template<int dim>
void test_renumbering(const FiniteElement<dim>& fe)
{
  std::vector<unsigned int> v(fe.dofs_per_cell);
  std::vector<std::vector<unsigned int> > start;
  FETools::compute_component_wise(fe, v, start);

  deallog << fe.get_name() << std::endl << '[';
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << v[i];
  deallog << " ]" << std::endl;

  for (unsigned int i=0;i<start.size();++i)
    {
      deallog << "Base " << i << ':';
      for (unsigned int j=0;j<start[i].size();++j)
	deallog << ' ' << start[i][j];
      deallog << std::endl;
    }
  
}


template<int dim>
void test_renumbering()
{
  deallog.push("Renumber");
  
  FE_Q<dim> q1(1);
  FE_Q<dim> q3(3);
//Todo:[GK] Test Raviart-Thomas and Nedelec?
  
  FESystem<dim> q1_3(q1,3);
  test_renumbering(q1_3);
  FESystem<dim> q3_3(q3,3);
  test_renumbering(q3_3);
  FESystem<dim> q3_2_q1_3(q3, 2, q1, 3);
  test_renumbering(q3_2_q1_3);

  deallog.pop();
}


int main()
{
  std::ofstream logfile("fe_tools.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test_projection<1>(logfile);
  test_projection<2>(logfile);
  test_projection<3>(logfile);
  
  test_renumbering<1>();
  test_renumbering<2>();
  test_renumbering<3>();
}
