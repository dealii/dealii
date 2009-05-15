/*
    Muti_Component_Function class implementation
*/
#include "../include/MultiComponentFunction.h"


// Constructor: Just call the Function constructor and set the component to 0
template<int dim> Multi_Component_Function<dim>::Multi_Component_Function( const double initial_time ):
                    Function<dim>( 1, initial_time ), component(0) {
}



// Set Component Function: Check that it is in range and then set it
template<int dim> void Multi_Component_Function<dim>::set_component(const unsigned int d ){
  // Check if the requested component is correct
  Assert( d<dim, ExcIndexRange( d, 0, dim ) );
  // We have only written 2d
  Assert( dim >= 2, ExcNotImplemented() );
  component = d;
}



// explicit template instantiation
template class Multi_Component_Function<deal_II_dimension>;
