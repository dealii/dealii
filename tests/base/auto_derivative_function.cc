//----------------------------  polynomial_test.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  polynomial_test.cc  ---------------------------

#include <iostream>


#include <cmath>

#include <base/point.h>
#include <base/auto_derivative_function.h>
#include <base/convergence_table.h>
#include <base/logstream.h>
#include <lac/vector.h>


template <int dim>
class AutoSinExp: public AutoDerivativeFunction<dim>
{
  public:
    AutoSinExp();
    virtual ~AutoSinExp() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

				     /**
				      * n_components=2. First
				      * component=0, second component
				      * filled by @p{value} function.
				      */
    virtual void vector_value (const Point<dim>   &p,
			       Vector<double>     &values) const;
};



template <int dim>
AutoSinExp<dim>::AutoSinExp():
		AutoDerivativeFunction<dim>(1e-6, 2)
{}



template <int dim>
double AutoSinExp<dim>::value (const Point<dim> &p,
			       const unsigned int) const
{
  return sin(2*p(0))*exp(3*p(1));
}


template <int dim>
void AutoSinExp<dim>::vector_value (const Point<dim> &p,
				    Vector<double>     &values) const
{
  Assert(values.size()==n_components, ExcDimensionMismatch(values.size(), n_components));
  values(0)=0;
  values(1)=value(p);
}

/*----------------------------------------------------------*/

template <int dim>
class ExactSinExp: public AutoSinExp<dim>
{
  public:
    ExactSinExp() {};
    ~ExactSinExp() {};
    
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				    const unsigned int  component = 0) const;    

    virtual void          vector_gradient (const Point<dim>            &p,
					   typename std::vector<Tensor<1,dim> > &gradients) const;
};




template <int dim>
Tensor<1,dim> ExactSinExp<dim>::gradient (const Point<dim> &p,
					 const unsigned int) const
{
  Tensor<1,dim> grad;
  grad[0]=2*cos(2*p(0))*exp(3*p(1));
  grad[1]=3*sin(2*p(0))*exp(3*p(1));
  return grad;
}


template <int dim>
void ExactSinExp<dim>::vector_gradient (const Point<dim>       &p,
				       typename std::vector<Tensor<1,dim> > &gradients) const
{
  Assert(gradients.size()==n_components, ExcDimensionMismatch(gradients.size(), n_components));
  
  gradients[0].clear();
  gradients[1]=gradient(p);
}


int main(int, char)
{
  std::ofstream logfile("auto_derivative_function.output");
  logfile.precision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);

  const unsigned int dim=2;
  AutoSinExp<dim> auto_function;
  ExactSinExp<dim> exact_function;
  Point<dim> p(0.23, 0.1);
  std::vector<Point<dim> > ps(1, p);
  
  Tensor<1,dim> u_grad=exact_function.gradient(p);


  AutoDerivativeFunction<dim>::DifferenceFormula formula;
  const double h_base=0.1;
  for (unsigned int order=1; order<5; ++order)
    {
      formula=AutoDerivativeFunction<dim>::get_formula_of_order(order);
      auto_function.set_formula(formula);
      deallog << "order=" << order << ",  formula=" << formula << std::endl;
      ConvergenceTable history;
      
      unsigned int factor=1;
      for (unsigned int i=0; i<6; ++i, factor*=2)
	{
	  history.add_value("f", factor);
	  history.omit_column_from_convergence_rate_evaluation("f");

	  auto_function.set_h(h_base/factor);

					   // Test of gradient function
	  Tensor<1,dim> a_grad=auto_function.gradient(p);
	  a_grad-=u_grad;
	  double value=sqrt(a_grad*a_grad);
	  history.add_value("grad", value);
	  history.set_scientific("grad", true);
	  history.set_precision("grad", 2);

					   // Test of gradient_list
					   // function
	  std::vector<Tensor<1,dim> > a_grads(1);
	  auto_function.gradient_list(ps, a_grads);
	  a_grads[0]-=u_grad;
	  value=sqrt(a_grads[0]*a_grads[0]);
	  history.add_value("grads[0]", value);
	  history.set_scientific("grads[0]", true);
	  history.set_precision("grads[0]", 2);
	  
					   // Test of vector_gradient
					   // function
	  std::vector<Tensor<1,dim> > a_vgrad(2);
	  auto_function.vector_gradient(p, a_vgrad);
	  a_vgrad[1]-=u_grad;
	  value=sqrt(a_vgrad[1]*a_vgrad[1]);
	  history.add_value("vgrad[1]", value);
	  history.set_scientific("vgrad[1]", true);
	  history.set_precision("vgrad[1]", 2);
	  
					   // Test of
					   // vector_gradient_list
					   // function
	  std::vector<std::vector<Tensor<1,dim> > >
	    a_vgrads(1, std::vector<Tensor<1,dim> > (2));
	  auto_function.vector_gradient_list(ps, a_vgrads);
	  a_vgrads[0][1]-=u_grad;
	  value=sqrt(a_vgrads[0][1]*a_vgrads[0][1]);
	  history.add_value("vgrads[1]", value);
	  history.set_scientific("vgrads[1]", true);
	  history.set_precision("vgrads[1]", 2);
	}
      history.evaluate_all_convergence_rates(
	ConvergenceTable::reduction_rate);
      history.write_text(deallog.get_file_stream());
    }
}
