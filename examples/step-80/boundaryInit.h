
//#----------------------------------------------------------
//#
//# This file defines the boundary and initial conditions
//#
//#----------------------------------------------------------

#ifndef GLOBAL_PARA
#define GLOBAL_PARA
#include "./globalPara.h"
#endif

//#----------------------------------------------------------
//# Declaration
//
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
class InitialValues : public Function<dim>
{
public:
  InitialValues () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};


//#----------------------------------------------------------
//# Implementation
//
template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &/*p*/,
                                   const unsigned int /*component*/) const
{
  return 293;
}


template <int dim>
double InitialValues<dim>::value (const Point<dim> &/*p*/,
                                   const unsigned int /*component*/) const
{
  return 293;
}





//#----------------------------------------------------------
//# Declaration and Implementation
//# mass density and heat capacity
//
template <int dim>
class RhoC : public Function<dim>
{
public:
  RhoC () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
double RhoC<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
//# p stores the xyz coordinates at each vertex
//# for 2D problems in xy, we assume the non-uniform is in y-axis.
  if ( p[1] >= -global_film_thickness )
  {
      return global_rho_Tio2 * global_C_Tio2;
  }
  else
      return global_rho_glass * global_C_glass;
}




//#----------------------------------------------------------
//# Declaration and Implementation
//# thermal conductivity
//
template <int dim>
class K_T : public Function<dim>
{
public:
  K_T () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;
};

template <int dim>
double K_T<dim>::value (const Point<dim> &p,
                                   const unsigned int /*component*/) const
{
//# p stores the xyz coordinates at each vertex
//# for 2D problems in xy, we assume the non-uniform is in y-axis.
  if ( p[1] >= -global_film_thickness)
  {
      return global_k_Tio2;
  }
  else
      return global_k_glass;
}

