	// C++ 
#include <iostream>

// The first base class:
#include "BaseClass.h"

// And now we define a second one:
template<int dim>
class Base2
{
public:
	Base2();
	~Base2();
	virtual int return_int(){return 321;};
};

template<int dim>
class Derived2 : public Base2<dim>
{
public:
	Derived2();
	~Derived2();
	virtual int return_int(){ return 123;}
};

template<int dim>
class DerivedDerived2 : public Derived2<dim>
{
public:
	DerivedDerived2();
	~DerivedDerived2();
	virtual int return_int(){ return 456;}
};

template<int dim>
Base2<dim>::Base2()
{ }

template<int dim>
Derived2<dim>::Derived2()
{ }

template<int dim>
DerivedDerived2<dim>::DerivedDerived2()
{ }

template<int dim>
Base2<dim>::~Base2()
{ }

template<int dim>
Derived2<dim>::~Derived2()
{ }

template<int dim>
DerivedDerived2<dim>::~DerivedDerived2()
{ }

/***************************************
 *
 ***************************************/
int main ( ) 
{
	const unsigned int dim = 2;
	
	//First test the dynamic library:
	{
	DerivedDerived<dim> *der = new DerivedDerived<dim>();
	
	const Base<dim>* my_class_base_pointer = der;
	
	if(dynamic_cast<const Derived<dim> *>(my_class_base_pointer) != 0)
		std::cout<<"SUCCESS"<<std::endl;
	else 
		std::cout<<"FAILURE"<<std::endl;	
	}
	
	//Now test the same exact code that was copied into this file:
	{
	DerivedDerived2<dim> *der = new DerivedDerived2<dim>();
	
	const Base2<dim>* my_class_base_pointer = der;
	
	if(dynamic_cast<const Derived2<dim> *>(my_class_base_pointer) != 0)
		std::cout<<"SUCCESS"<<std::endl;
	else 
		std::cout<<"FAILURE"<<std::endl;	
	}
  
  return 0;
}


