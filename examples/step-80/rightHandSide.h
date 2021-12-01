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
class RightHandside : public Function<dim>
{
public:
  RightHandside () : Function<dim>() {}

  double value_v2 (const Point<dim>   &p);
};

//#----------------------------------------------------------
//# Implementation
//
template <int dim>
double RightHandside<dim>::value_v2 (const Point<dim> &p)
{

    double Rnp_inter;

    double alpha_abs = 1e4;


    if(p[1] >= -global_film_thickness)
    {

    double P00 = global_Pow_laser / global_PI / global_c_laser / global_c_laser / 2.0;

    double I00 = P00 * std::exp(-
                ( 
                    (p[0] - global_V_scan_x * this->get_time()-global_init_position_x0) * 
                    (p[0] - global_V_scan_x * this->get_time()-global_init_position_x0)
                ) /
                (2.0 * global_c_laser * global_c_laser) );


    return alpha_abs * I00 *  std::exp(-alpha_abs*(0 - p[1]));

    }

    else
    {
        //# no absorption
        return 0;
    }

}


