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
#include <base/polynomial.h>



bool equals_delta_ij(double value, unsigned int i, unsigned int j)
{
  double eps=1e-14;
  if ((i==j && fabs(value-1)<eps) || (i!=j && fabs(value)<eps))
    return true;
  else
    return false;
}




int main(int, char)
{
  ofstream logfile("polynomial_test.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  vector<double> values(1);
  deallog << "LagrangeEquidistant polynoms:" << endl;
  for (unsigned int order=1; order<=4; ++order)
    {
      deallog << "Polynomial p of order " << order << endl;
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
	      deallog << endl;

					       // now also check
					       // whether the other
					       // @p{value} function
					       // returns the same
					       // result
	      if (polynom.value(x) != values[0])
		{
		  deallog << "The two `value' functions return different results!"
			  << endl;
		  abort ();
		};
	    }
	}
    }

  deallog << endl << "Test derivatives computed by the Horner scheme:" << endl;
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
//    	       << "   v_exact[i]=" << v_exact[i] << endl;
	  if (fabs(v_horner[i]-v_exact[i])>1e-12)
	    ok=false;
	}

      if (ok)
	deallog << "ok";
      else
	deallog << "false";

      deallog << endl;
    }
}



