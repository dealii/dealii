#include <base/tensor.h>
#include <iostream>

int main () {
  const unsigned int dim=3;
  
  Tensor<2,dim> t;
  Tensor<2,dim> tt;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      cout << i << ' ' << j << ':' << t[i][j] << endl;
  cout << "----------------------------" << endl;
	
  t[0][0] = 1;
  t[0][1] = 2;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      cout << i << ' ' << j << ':' << t[i][j] << endl;
  cout << "----------------------------" << endl;

  contract (tt,t,t);

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      cout << i << ' ' << j << ':' << tt[i][j] << endl;
};
