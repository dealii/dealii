/*---------------------------- functiontime.h ---------------------------*/
/*      $Id$ */
/*      Copyright W. Bangerth, Guido Kanschat
 * University of Heidelberg, 1998, 1999 */
#ifndef __functiontime_H
#define __functiontime_H
/*---------------------------- functiontime.h ---------------------------*/


#include <base/exceptions.h>

/**
 *  Support for time dependent functions.
 *  The library was also designed for time dependent problems. For this
 *  purpose, the function objects also contain a field which stores the
 *  time, as well as functions manipulating them. Time independent problems
 *  should not access or even abuse them for other purposes, but since one
 *  normally does not create thousands of function objects, the gain in
 *  generality weighs out the fact that we need not store the time value
 *  for not time dependent problems. The second advantage is that the derived
 *  standard classes like #ZeroFunction#, #ConstantFunction# etc also work
 *  for time dependent problems.
 *
 *  Access to the time goes through the following functions:
 *  \begin{verbatim}
 *  \item #get_time#: return the present value of the time variable.
 *  \item #set_time#: set the time value to a specific value.
 *  \item #advance_time#: increase the time by a certain time step.
 *  \end{verbatim}
 *  The latter two functions are virtual, so that derived classes can
 *  perform computations which need only be done once for every new time.
 *  For example, if a time dependent function had a factor #sin(t)#, then
 *  it may be a reasonable choice to calculate this factor in a derived
 *  version of #set_time#, store it in a member variable and use that one
 *  rather than computing it every time #operator()#, #value_list# or one
 *  of the other functions is called.
 *
 *  By default, the #advance_time# function calls the #set_time# function
 *  with the new time, so it is sufficient in most cases to overload only
 *  #set_time# for computations as sketched out above.
 *
 *  The constructor of this class takes an initial value for the time
 *  variable, which defaults to zero. Because a default value is given,
 *  none of the derived classes needs to take an initial value for the
 *  time variable if not needed.
 *
 *  Once again the warning: do not use the #time# variable for any other
 *  purpose than the intended one! This will inevitably lead to confusion.
 *
 *
 *  @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999
 */
class FunctionTime {
  public:
				     /**
				      * Constructor. May take an initial vakue
				      * for the time variable, which defaults
				      * to zero.
				      */
    FunctionTime (const double initial_time = 0.0);

				     /**
				      * Virtual destructor.
				      */
    virtual ~FunctionTime();
  
				     /**
				      * Return the value of the time variable/
				      */
    double get_time () const;

				     /**
				      * Set the time to #new_time#, overwriting
				      * the old value.
				      */
    virtual void set_time (const double new_time);

				     /**
				      * Advance the time by the given
				      * time step #delta_t#.
				      */
    virtual void advance_time (const double delta_t);

  private:
				     /**
				      * Store the present time.
				      */
    double time;
};

/*------------------------------ Inline functions ------------------------------*/

inline double
FunctionTime::get_time () const
{
  return time;
};

/*---------------------------- functiontime.h ---------------------------*/
/* end of #ifndef __functiontime_H */
#endif
/*---------------------------- functiontime.h ---------------------------*/
