#include "BaseClass.h"

template<int dim>
Base<dim>::Base()
{ }

template<int dim>
Derived<dim>::Derived()
{ }

template<int dim>
DerivedDerived<dim>::DerivedDerived()
{ }

template<int dim>
Base<dim>::~Base()
{ }

template<int dim>
Derived<dim>::~Derived()
{ }

template<int dim>
DerivedDerived<dim>::~DerivedDerived()
{ }

template class DerivedDerived<3>;
template class DerivedDerived<2>;
template class DerivedDerived<1>;
template class DerivedDerived<0>;

template class Derived<3>;
template class Derived<2>;
template class Derived<1>;
template class Derived<0>;

template class Base<3>;
template class Base<2>;
template class Base<1>;
template class Base<0>;