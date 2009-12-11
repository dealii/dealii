#include "BaseClass.h"

int main ( ) 
{
  const unsigned int dim = 2;
	
  DerivedDerived<dim> *der = new DerivedDerived<dim>();
  
  const Base<dim>* my_class_base_pointer = der;
  
  if(dynamic_cast<const Derived<dim> *>(my_class_base_pointer) != 0)
    return 0;
  else 
    return 1;
}


