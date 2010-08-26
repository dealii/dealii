//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// Check methods of the Q1WedgeFunction
// 1-dim does not exist
// 2-dim Q1WedgeFunction = xy
// 3-dim Q1WedgeFunction = xy

#include "../tests.h"
#include <base/data_out_base.h>
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <base/function_lib.h>
#include <base/auto_derivative_function.h>
#include <base/flow_function.h>
#include <lac/vector.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>


template <int dim>
void check_values(const Function<dim>& f)
{
  Point<dim> p;
  for(unsigned int i=0; i<dim; ++i)
    p[i] = i+1;
  deallog << f.value(p) << std::endl; 
  deallog << " values checked" << std::endl;
}

template <int dim>
void check_value_list(const Function<dim>& f)
{
  const unsigned int max_points = 3;
  std::vector< Point<dim> >  point_vector(max_points);
  for(unsigned int j=0; j<max_points; ++j)
  {
    Point<dim> p;
    for(unsigned int i=0; i<dim; ++i)
      p[i]=j+1;
    point_vector[j] = p;
  }
  std::vector< double > values(max_points);
  f.value_list(point_vector, values);
  for(unsigned int j=0; j<max_points; ++j)
    deallog << values[j] << std::endl; 
  deallog << " value_list checked" << std::endl;
}

template <int dim>
void check_vector_value_list(const Function<dim>& f)
{
  const unsigned int max_points = 3;
  std::vector< Point<dim> >  point_vector(max_points);
  for(unsigned int j=0; j<max_points; ++j)
  {
    Point<dim> p;
    for(unsigned int i=0; i<dim; ++i)
      p[i]=j+1;
    point_vector[j] = p;
  }
  std::vector< Vector<double> > values(max_points, Vector<double>(1));
  f.vector_value_list(point_vector, values);
  for(unsigned int j=0; j<max_points; ++j)
    deallog << values[j](0) << std::endl; 
  deallog << " vector_value_list checked" << std::endl;
}

template <int dim>
void check_gradients(const Function<dim>& f)
{
  Point<dim> p;
  for(unsigned int i=0; i<dim; ++i)
  {
    p[i] = i+3;
  }
  Tensor <1, dim> grads = f.gradient(p);
  for(unsigned int i=0; i<dim; ++i)
  {
    deallog << i <<"-der: " << grads[i] << std::endl; 
  }
  deallog << " gradients checked" << std::endl;
}

template <int dim>
void check_gradient_list(const Function<dim>& f)
{
  const unsigned int max_points = 3;
  std::vector< Point<dim> >  point_vector(max_points);
  for(unsigned int j=0; j<max_points; ++j)
  {
    Point<dim> p;
    for(unsigned int i=0; i<dim; ++i)
      p[i]=j+1;
    point_vector[j] = p;
  }
  std::vector< Tensor<1,dim> > grads(max_points);
  f.gradient_list(point_vector, grads);
  for(unsigned int j=0; j<max_points; ++j)
    deallog << grads[j] << std::endl; 
  deallog << " gradient_list checked" << std::endl;
}

template <int dim>
void check_vector_gradient_list(const Function<dim>& f)
{
  const unsigned int max_points = 3;
  std::vector< Point<dim> >  point_vector(max_points);
  for(unsigned int j=0; j<max_points; ++j)
  {
    Point<dim> p;
    for(unsigned int i=0; i<dim; ++i)
      p[i]=j+1;
    point_vector[j] = p;
  }
  std::vector< std::vector< Tensor<1,dim> > >  gradients(max_points, 
            std::vector< Tensor<1,dim> > (1));
  f.vector_gradient_list(point_vector, gradients);
  for(unsigned int j=0; j<max_points; ++j)
    deallog << gradients[j][0] << std::endl; 
  deallog << " vector_gradient_list checked" << std::endl;
}

template <int dim>
void check_laplacian(const Function<dim>& f)
{
  Point<dim> p;
  for(unsigned int i=0; i<dim; ++i)
    p[i] = i+1;
  Assert( std::fabs(f.laplacian(p)) < 1e-12,
         ExcInternalError());
}
  
template <int dim>
void check_laplacian_list(const Function<dim>& f)
{
  const unsigned int max_points = 3;
  std::vector< Point<dim> >  point_vector(max_points);
  for(unsigned int j=0; j<max_points; ++j)
  {
    Point<dim> p;
    for(unsigned int i=0; i<dim; ++i)
      p[i]=j+1;
    point_vector[j] = p;
  }
  std::vector< double > values(max_points);
  f.laplacian_list(point_vector, values);
  for(unsigned int j=0; j<max_points; ++j)
    Assert( std::fabs(values[j]) < 1e-12,
         ExcInternalError());
}

int main()
{
  std::ofstream logfile("functions_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  if (true)
    {
      deallog << " Functions::Q1WedgeFunction<2>" << std::endl;
      check_values(Functions::Q1WedgeFunction<2>());
      check_values(Functions::Q1WedgeFunction<3>());
      check_value_list(Functions::Q1WedgeFunction<2>());
      check_value_list(Functions::Q1WedgeFunction<3>());
      check_vector_value_list(Functions::Q1WedgeFunction<2>());
      check_vector_value_list(Functions::Q1WedgeFunction<3>());
      check_gradients(Functions::Q1WedgeFunction<2>());
      check_gradients(Functions::Q1WedgeFunction<3>());
      check_gradient_list(Functions::Q1WedgeFunction<2>());
      check_gradient_list(Functions::Q1WedgeFunction<3>());
      check_vector_gradient_list(Functions::Q1WedgeFunction<2>());
      check_vector_gradient_list(Functions::Q1WedgeFunction<3>());
      check_laplacian(Functions::Q1WedgeFunction<2>());
      check_laplacian(Functions::Q1WedgeFunction<3>());
    }

}
