/*----------------------------   function.h     ---------------------------*/
/*      $Id$                 */
#ifndef __function_H
#define __function_H
/*----------------------------   function.h     ---------------------------*/


//forward declaration
template <int dim> class Point;
template <class T, class Alloc=alloc> class vector;




template <int dim>
class Function {
  public:
    virtual double operator () (const Point<dim> &p) const = 0;
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values) const;
};




/*----------------------------   function.h     ---------------------------*/
/* end of #ifndef __function_H */
#endif
/*----------------------------   function.h     ---------------------------*/
