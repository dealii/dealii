/*
    Header file for the classes that define the exact solution
    and driving force.

    These classes are much like any other Function class
    so we are not going to comment on them

    by Abner Salgado.
*/
#ifndef _EQ_DATA_H_
#define _EQ_DATA_H_


/*
    We want to be able to select which component of vector
    valued functions we are going to work on
*/
#include "../include/MultiComponentFunction.h"


// this is dealii
#include <base/point.h>


// We need to define sines and cosines
#include <cmath>



// This is ugly, we have to remove it
const double PI = std::acos( -1. );



// The Velocity function
template<int dim> class Velocity: public Multi_Component_Function<dim> {
  public:
    Velocity(const double initial_time =0.0);
    virtual double value(const Point<dim> &p, const unsigned int component = 0) const;
    virtual Tensor<1,dim> gradient(const Point<dim> &p, const unsigned int component=0) const;
    virtual void value_list( const std::vector< Point<dim> > &points, std::vector<double> &values,
                               const unsigned int component = 0 ) const;
    virtual void gradient_list( const std::vector< Point<dim> > &points, std::vector< Tensor<1,dim> > &gradients,
                                  const unsigned int component = 0 ) const;
};



// The Pressure function
template<int dim> class Pressure: public Function<dim>{
  public:
    Pressure(const double initial_time = 0.0);
    virtual double value(const Point<dim> &p, const unsigned int component = 0) const;
    virtual Tensor<1,dim> gradient(const Point<dim> &p, const unsigned int component=0) const;
    virtual void value_list( const std::vector< Point<dim> > &points, std::vector<double> &values,
                              const unsigned int component = 0 ) const;
    virtual void gradient_list( const std::vector< Point<dim> > &points, std::vector< Tensor<1,dim> > &gradients,
                                  const unsigned int component = 0 ) const;
};



// The Force function
template<int dim> class Force: public Multi_Component_Function<dim>{
  public:
    Force( const double initial_time =0.0 );
    virtual double value( const Point<dim> &p, const unsigned int component = 0 ) const;
    virtual void value_list( const std::vector< Point<dim> > &points, std::vector<double> &values, const unsigned int component = 0 ) const;
};

#endif
