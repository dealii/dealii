/*----------------------------   fe_lib.h     ---------------------------*/
/*      $Id$                 */
#ifndef __fe_lib_H
#define __fe_lib_H
/*----------------------------   fe_lib.h     ---------------------------*/


#include <fe/fe.h>



/**
  Define a (bi-, tri-, etc)linear finite element in #dim# space dimensions.
  */
template <int dim>
class FELinear : public FiniteElement<dim> {
  public:
    FELinear ();
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;
};




/**
  Define a (bi-, tri-, etc)quadratic finite element in #dim# space dimensions.
  */
template <int dim>
class FEQuadratic : public FiniteElement<dim> {
  public:
    FEQuadratic ();
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;
};




/**
  Define a (bi-, tri-, etc)cubic finite element in #dim# space dimensions.
  */
template <int dim>
class FECubic : public FiniteElement<dim> {
  public:
    FECubic ();
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;
};



/*----------------------------   fe_lib.h     ---------------------------*/
/* end of #ifndef __fe_lib_H */
#endif
/*----------------------------   fe_lib.h     ---------------------------*/
