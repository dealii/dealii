// $Id$

#include <base/tensor.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <fstream>


int main ()
{
  ofstream logfile("tensor.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  double a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};
  double b[3][3] = {{25,31,37}, {45,57,69}, {75,96,117}};
    
  const unsigned int dim=3;
  Tensor<2,dim> t(a);
  Tensor<2,dim> tt;
  Tensor<2,dim> result(b);

  Vector<double> unrolled(9);
  
  t.unroll(unrolled);
  deallog << "unrolled:";
  for (unsigned i=0;i<9;i++)
    deallog << ' ' << unrolled(i);
  deallog << endl;

  deallog << "t=" << endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
	deallog << t[i][j] << ' ';
      deallog << endl;
    };
  deallog << endl;
  
  contract (tt,t,t);

  deallog << "tt=" << endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
	deallog << tt[i][j] << ' ';
      deallog << endl;
    };
  deallog << endl;

  if (tt==result)
    {
      deallog << "Result OK." << endl;
      return 0;
    }
  else
    {
      deallog << "Result WRONG!" << endl;
      return 1;
    };
};
