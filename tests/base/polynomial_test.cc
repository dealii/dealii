//----------------------------  polynomial_test.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  polynomial_test.cc  ---------------------------

#include <iostream>
#include <fstream>

#include <base/logstream.h>
#include <base/tensor_product_polynomials.h>



bool equals_delta_ij(double value, unsigned int i, unsigned int j)
{
  double eps=1e-14;
  if ((i==j && fabs(value-1)<eps) || (i!=j && fabs(value)<eps))
    return true;
  else
    return false;
}

void Q3_4th_shape_function_values_and_grads_dim2(
  const Point<2> &point, double &v_exact,
  Tensor<1,2> &grad_exact, Tensor<2,2> &grad2);




int main(int, char)
{
  ofstream logfile("polynomial_test.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  vector<double> values(1);
  deallog << "LagrangeEquidistant polynoms:" << std::endl;
  for (unsigned int order=1; order<=4; ++order)
    {
      deallog << "Polynomial p of order " << order << std::endl;
      for (unsigned int s_point=0; s_point<=order; ++s_point)
	{
	  LagrangeEquidistant polynom(order, s_point);

					   // support points in vertices
	  for (unsigned int i=0; i<=order; ++i)
	    {
	      double x=static_cast<double>(i)/order;
	      polynom.value(x, values);
	      deallog << " p_" << s_point << "(" << x << ")";
//	      deallog << "=" << values[0];
	      if (equals_delta_ij(values[0], s_point, i))
		deallog << "   ok";
	      else
		deallog << "   false";
	      deallog << std::endl;

					       // now also check
					       // whether the other
					       // @p{value} function
					       // returns the same
					       // result
	      if (polynom.value(x) != values[0])
		{
		  deallog << "The two `value' functions return different results!"
			  << std::endl;
		  abort ();
		};
	    }
	}
    }

  deallog << std::endl << "Test derivatives computed by the Horner scheme:" << std::endl;
  LagrangeEquidistant pol(4, 2);
  vector<double> v_horner(6);
  for (unsigned int i=0; i<=10; ++i)
    {
      double xi=i*0.1;
      deallog << "x=" << xi << ",    all derivatives: ";
      vector<double> v_exact(6);
      
      v_exact[1]=256.0*xi*xi*xi-384.0*xi*xi+152.0*xi-12.0;
      v_exact[0]=64.0*xi*xi*xi*xi-128.0*xi*xi*xi+76.0*xi*xi-12.0*xi;
      v_exact[2]=768.0*xi*xi-768.0*xi+152.0;
      v_exact[3]=1536*xi-768;
      v_exact[4]=1536;
      v_exact[5]=0;

      pol.value(xi, v_horner);

      bool ok=true;
      for (unsigned int i=0; i<v_exact.size(); ++i)
	{
//    	  deallog << "v_horner[i]=" << v_horner[i]
//    	       << "   v_exact[i]=" << v_exact[i] << std::endl;
	  if (fabs(v_horner[i]-v_exact[i])>1e-12)
	    ok=false;
	}

      if (ok)
	deallog << "ok";
      else
	deallog << "false";

      deallog << std::endl;
    }


  deallog << std::endl << "Test of TensorProductPolynomials:" << std::endl;
  deallog << "2D Example:" << std::endl;
  unsigned int p=3,
   n_tensor_pols=(p+1)*(p+1);
  vector<SmartPointer<Polynomial> > pols;
  
  for (unsigned int i=0; i<=p; ++i)
    pols.push_back(new LagrangeEquidistant(p, i));
  
  TensorProductPolynomials<2> tp_pol(pols);
  vector<double> vs(n_tensor_pols);
  vector<Tensor<1,2> > grads(n_tensor_pols);
  vector<Tensor<2,2> > grad_grads(n_tensor_pols);

  double v_exact;
  Tensor<1,2> grad_exact;
  Tensor<2,2> grad_grad_exact;
  
  double xi=0.35;
  double eta=0.62;
  Point<2> point(xi,eta);
  tp_pol.compute(point, vs, grads, grad_grads);

				   // 4th shape function of Q3<2> is
				   // equivalent to its 1st shape
				   // function in lexicographical
				   // order.
  Q3_4th_shape_function_values_and_grads_dim2(point, v_exact, grad_exact, grad_grad_exact);

  unsigned int i=1;
  deallog << "v_" << i << "=" << vs[i] << std::endl;
  deallog << "v_exact=" << v_exact << std::endl;
  deallog << "grad_v_" << i << "=" << grads[i] << std::endl;
  deallog << "grad_exact=" << grad_exact << std::endl;
  for (unsigned int j=0; j<grad_grads[i].dimension; ++j)
    for (unsigned int k=0; k<grad_grads[i].dimension; ++k)
      {
	deallog << "grad_grads_" << i<< "[" << j << "][" << k << "]="
		<< grad_grads[i][j][k] << std::endl;
	deallog << "grad2_exact[" << j << "][" << k << "]="
		<< grad_grad_exact[j][k] << std::endl;
      }
}



void Q3_4th_shape_function_values_and_grads_dim2(
  const Point<2> &point, double &v_exact, Tensor<1,2> &grad_exact, Tensor<2,2> &grad2)
{
				   // the following functions
				   // are taken from fe_lib.cubic.cc
  const double xi=point(0),
	      eta=point(1);
  
  v_exact=9.0*xi-45.0/2.0*xi*xi+27.0/2.0*xi*xi*xi+(
    -99.0/2.0*xi+495.0/4.0*xi*xi-297.0/4.0*xi*xi*xi)*eta+(81.0*xi-405.0/2.0*xi*xi+243.0/2.0*xi*xi*xi)*
		 eta*eta+(-81.0/2.0*xi+405.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta*eta*eta;
  
  grad_exact[0]=9.0-45.0*xi+81.0/2.0*xi*xi+(-99.0/2.0+495.0/2.0*xi-891.0/4.0*xi*xi)*eta+
		(81.0-405.0*xi+729.0/2.0*xi*xi)*eta*eta+(-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta*eta;
  grad_exact[1]=-99.0/2.0*xi+495.0/4.0*xi*xi-297.0/4.0*xi*xi*xi+2.0*(
    81.0*xi-405.0/2.0*xi*xi+243.0/2.0*xi*xi*xi)*eta+3.0*(
      -81.0/2.0*xi+405.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta*eta;
  
  grad2[0][0] = -45.0+81.0*xi+(495.0/2.0-891.0/2.0*xi)*eta+(-405.0+729.0*xi)*eta*eta+
		(405.0/2.0-729.0/2.0*xi)*eta*eta*eta;
  grad2[0][1] = -99.0/2.0+495.0/2.0*xi-891.0/4.0*xi*xi+2.0*(81.0-405.0*xi+729.0/2.0*xi*xi)*eta+
		3.0*(-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta;
  grad2[1][0] = -99.0/2.0+495.0/2.0*xi-891.0/4.0*xi*xi+2.0*(81.0-405.0*xi+729.0/2.0*xi*xi)*eta+
		3.0*(-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta;
  grad2[1][1] = 162.0*xi-405.0*xi*xi+243.0*xi*xi*xi+6.0*(-81.0/2.0*xi+405.0/4.0*xi*xi-
							 243.0/4.0*xi*xi*xi)*eta;
}

