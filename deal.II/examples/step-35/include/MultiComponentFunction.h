/*
    Header file for Multi Component Functions.
    The whole purpose of this class to exist is so we do not have to write the function
        void set_component( const unsigned int d )
    twice

    by Abner Salgado.
*/
#ifndef _MULTI_COMPONENT_FUNCTION_H_
#define _MULTI_COMPONENT_FUNCTION_H_


/*
  This is basically a workaround for vector valued functions
  when we just want one of its components
  so they are derived from the Function class
*/
#include <base/function.h>


using namespace dealii;


/*
    The only funcitonality this class has is that it provides a common wrapper
    for vector valued functions for the case when you only want to ask them for
    one of their components
*/
template<int dim> class Multi_Component_Function: public Function<dim> {
  public:
    /*
      Constructor. It does not do anything interesting but set the initial
      time for the function to evaluate its values
    */
    Multi_Component_Function( const double initial_time = 0. );
    /*
      This is the whole reason this class exists.
      It wouldn't seem logical to have this written twice. One for the
      Velocity function and one for the Force function.
      So we are instead going to derive these from this class
      and so adding this functionality of selecting which component you
      want to give the value of
    */
    void set_component(const unsigned int d );
  protected:
    // The current component we are working on
    unsigned int component;
};


#endif
