#include <base/tensor.h>
#include <iostream>


DeclException2 (Exc1,
		int, int,
		<< "arg1=" << arg1 << " arg2=" << arg2);


int main () {
  double a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};
  double b[3][3] = {{25,31,37}, {45,57,69}, {75,96,117}};

  try 
    {
      ThrowIfNot (1>2, Exc1(1,2));
      cerr << "11111111111111111111";
    }
  catch (exception &e)
    {
      cerr << e.what();
      cerr << "33333333333333333333";
    };
  
    
  const unsigned int dim=3;
  Tensor<2,dim> t(a);
  Tensor<2,dim> tt;
  Tensor<2,dim> result(b);

  cout << "t=" << endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
	cout << t[i][j] << ' ';
      cout << endl;
    };
  cout << endl;
  
  contract (tt,t,t);

  cout << "tt=" << endl;
  for (unsigned int i=0; i<dim; ++i)
    {
      for (unsigned int j=0; j<dim; ++j)
	cout << tt[i][j] << ' ';
      cout << endl;
    };
  cout << endl;

  if (tt==result)
    {
      cout << "Result OK." << endl;
      return 0;
    }
  else
    {
      cout << "Result WRONG!" << endl;
      return 1;
    };
};
