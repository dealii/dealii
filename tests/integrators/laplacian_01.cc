//--------------------------------------------------------------------
//    $Id: cochain_01.cc 25686 2012-07-04 14:18:43Z bangerth $
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------

// Test the functions in integrators/laplace.h
// Output matrices and assert consistency of residuals

#include "../tests.h"
#include "../lib/test_grids.h"

#include <deal.II/base/logstream.h>
#include <deal.II/integrators/laplace.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>

using namespace LocalIntegrators::Laplace;

template <int dim>
void test_cell(const FEValuesBase<dim>& fev)
{
  const unsigned int n = fev.dofs_per_cell;
  unsigned int d=fev.get_fe().n_components();
  FullMatrix<double> M(n,n);
  cell_matrix(M,fev);
  M.print(deallog,8);
  
  Vector<double> u(n), v(n), w(n);
  std::vector<std::vector<Tensor<1,dim> > >
    ugrad(d,std::vector<Tensor<1,dim> >(fev.n_quadrature_points));
  
  std::vector<unsigned int> indices(n);
  for (unsigned int i=0;i<n;++i)
    indices[i] = i;
  
  deallog << "Residuals";
  for (unsigned int i=0;i<n;++i)
    {
      u = 0.;
      u(i) = 1.;
      w = 0.;
      fev.get_function_gradients(u, indices, ugrad, true);
      cell_residual(w, fev, make_slice(ugrad));
      M.vmult(v,u);
      w.add(-1., v);
      deallog << ' ' << w.l2_norm();
      if (d==1)
	{
	  cell_residual(w, fev, ugrad[0]);
	  M.vmult(v,u);
	  w.add(-1., v);
	  deallog << " e" << w.l2_norm();	  
	}
    }
  deallog << std::endl;
}


template <int dim>
void test_boundary(const FEValuesBase<dim>& fev)
{
  const unsigned int n = fev.dofs_per_cell;
  unsigned int d=fev.get_fe().n_components();
  FullMatrix<double> M(n,n);
  nitsche_matrix(M, fev, 17);
  M.print(deallog,8);
  
  Vector<double> u(n), v(n), w(n);
  std::vector<std::vector<double> >
    uval    (d,std::vector<double>(fev.n_quadrature_points)),
    null_val(d,std::vector<double>(fev.n_quadrature_points, 0.));
  std::vector<std::vector<Tensor<1,dim> > >
    ugrad   (d,std::vector<Tensor<1,dim> >(fev.n_quadrature_points));
  
  std::vector<unsigned int> indices(n);
  for (unsigned int i=0;i<n;++i)
    indices[i] = i;
  
  deallog << "Residuals";
  for (unsigned int i=0;i<n;++i)
    {
      u = 0.;
      u(i) = 1.;
      w = 0.;
      fev.get_function_values(u, indices, uval, true);
      fev.get_function_gradients(u, indices, ugrad, true);
      nitsche_residual(w, fev, make_slice(uval), make_slice(ugrad), make_slice(null_val), 17);
      M.vmult(v,u);
      w.add(-1., v);
      deallog << ' ' << w.l2_norm();
      if (d==1)
	{
	  nitsche_residual(w, fev, uval[0], ugrad[0], null_val[0], 17);
	  M.vmult(v,u);
	  w.add(-1., v);
	  deallog << " e" << w.l2_norm();	  
	}
    }
  deallog << std::endl;
}


template <int dim>
void test_face(const FEValuesBase<dim>& fev1,
	       const FEValuesBase<dim>& fev2)
{
  const unsigned int n1 = fev1.dofs_per_cell;
  const unsigned int n2 = fev1.dofs_per_cell;
  FullMatrix<double> M11(n1,n1);
  FullMatrix<double> M12(n1,n2);
  FullMatrix<double> M21(n2,n1);
  FullMatrix<double> M22(n2,n2);
//  Vector<double> u(n), v(n);
  
  ip_matrix(M11, M12, M21, M22, fev1, fev2, 17);
  deallog << "M11" << std::endl;
  M11.print(deallog,8);
  deallog << "M22" << std::endl;
  M22.print(deallog,8);
  deallog << "M12" << std::endl;
  M12.print(deallog,8);
  deallog << "M21" << std::endl;
  M21.print(deallog,8);
}


template <int dim>
void
test_fe(Triangulation<dim>& tr, FiniteElement<dim>& fe)
{
  deallog << fe.get_name() << std::endl << "cell matrix" << std::endl;
  QGauss<dim> quadrature(fe.tensor_degree()+1);
  FEValues<dim> fev(fe, quadrature, update_gradients);

  typename Triangulation<dim>::cell_iterator cell1 = tr.begin(1);
  fev.reinit(cell1);
  test_cell(fev);
  
  QGauss<dim-1> face_quadrature(fe.tensor_degree()+1);
  FEFaceValues<dim> fef1(fe, face_quadrature, update_values | update_gradients | update_normal_vectors);
  for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
    {
      deallog << "boundary_matrix " << i << std::endl;
      fef1.reinit(cell1, i);
      test_boundary(fef1);
    }

  FEFaceValues<dim> fef2(fe, face_quadrature, update_values | update_gradients);
  typename Triangulation<dim>::cell_iterator cell2 = cell1->neighbor(1);
  
  deallog << "face_matrix " << 0 << std::endl;
  cell1 = tr.begin(1);
  fef1.reinit(cell1, 1);
  fef2.reinit(cell2, 0);
  test_face(fef1, fef2);
  
  
}


template <int dim>
void
test(Triangulation<dim>& tr)
{
  FE_DGQ<dim> q1(1);
  test_fe(tr, q1);

//  FE_DGQ<dim> q2(2);
//  test_fe(tr, q2);

 FE_Nedelec<dim> n1(1);  
 test_fe(tr, n1);  
}


int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  Triangulation<2> tr2;
  TestGrids::hypercube(tr2, 1);
  test(tr2);
  
//  Triangulation<3> tr3;
//  TestGrids::hypercube(tr3, 1);
//  test(tr3);
  
}
