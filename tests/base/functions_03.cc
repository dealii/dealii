//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// Checking the PillowFunction

#include "../tests.h"
#include <base/data_out_base.h>
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <base/function_lib.h>
#include <lac/vector.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>



template <int dim>
void check_value (const Function<dim>& f)
{
  Point<dim> p;
  for (unsigned int i=0; i<dim; i++)
    p[i] = i;
  
  deallog << f.value(p) << std::endl;
  deallog << "values checked" << std::endl;

}

template <int dim>
void check_value_list (const Function<dim>& f)
{
  const unsigned int max_number_of_points = 5;
  std::vector< Point<dim> >  points(max_number_of_points);
  
  for(unsigned int i=0; i<max_number_of_points; ++i)
  {
    Point<dim> p;
    for(unsigned int j=0; j<dim; ++j)
      p[j]=i+1;
    
    points[i] = p;
  }
  std::vector< double > values(max_number_of_points);
  f.value_list(points, values);

  for(unsigned int j=0; j<max_number_of_points; ++j)
    deallog << values[j] << std::endl; 

  deallog << " value_list checked" << std::endl;
}


template <int dim>
void check_gradient (const Function<dim>& f)
{
  Point<dim> p;
  for (unsigned int i=0; i<dim; i++)
    p[i] = i;
  
  deallog << f.gradient(p) << std::endl;
  deallog << " gradients checked" << std::endl;
}

template <int dim>
void check_gradient_list (const Function<dim>& f)
{
  const unsigned int max_number_of_points = 5;
  std::vector< Point<dim> >  points(max_number_of_points);
  
  for(unsigned int i=0; i<max_number_of_points; ++i)
  {
    Point<dim> p;
    for(unsigned int j=0; j<dim; ++j)
      p[j]=i+1;
    
    points[i] = p;
  }
  std::vector< Tensor<1,dim> > tensors(max_number_of_points);
  f.gradient_list(points, tensors);

  for(unsigned int j=0; j<max_number_of_points; ++j)
    deallog << tensors[j] << std::endl; 

  deallog << " gradient_list checked" << std::endl;
}


template <int dim>
void check_laplacian (const Function<dim>& f)
{
  Point<dim> p;
  for (unsigned int i=0; i<dim; i++)
    p[i] = i;
  
  deallog << f.laplacian(p) << std::endl;
  deallog << " laplacians checked" << std::endl;
}


template <int dim>
void check_laplacian_list (const Function<dim>& f)
{
  const unsigned int max_number_of_points = 5;
  std::vector< Point<dim> >  points(max_number_of_points);
  
  for(unsigned int i=0; i<max_number_of_points; ++i)
  {
    Point<dim> p;
    for(unsigned int j=0; j<dim; ++j)
      p[j]=i+1;
    
    points[i] = p;
  }
  std::vector< double > values(max_number_of_points);
  f.laplacian_list(points, values);

  for(unsigned int j=0; j<max_number_of_points; ++j)
    deallog << values[j] << std::endl; 
 
  deallog << " laplacian_list checked" << std::endl;
}


int main()
{
  std::ofstream logfile("functions_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  deallog << "Functions PillowFunction" << std::endl;
  check_value (Functions::PillowFunction<1>());
  check_value (Functions::PillowFunction<2>());
  check_value (Functions::PillowFunction<3>());
  
  check_value_list (Functions::PillowFunction<1>());
  check_value_list (Functions::PillowFunction<2>());
  check_value_list (Functions::PillowFunction<3>());
  
  check_gradient (Functions::PillowFunction<1>());
  check_gradient (Functions::PillowFunction<2>());
  check_gradient (Functions::PillowFunction<3>());
  
  check_gradient_list (Functions::PillowFunction<1>());
  check_gradient_list (Functions::PillowFunction<2>());
  check_gradient_list (Functions::PillowFunction<3>());
  
  check_laplacian (Functions::PillowFunction<1>());
  check_laplacian (Functions::PillowFunction<2>());
  check_laplacian (Functions::PillowFunction<3>());
  
  check_laplacian_list (Functions::PillowFunction<1>());
  check_laplacian_list (Functions::PillowFunction<2>());
  check_laplacian_list (Functions::PillowFunction<3>());

     
}



