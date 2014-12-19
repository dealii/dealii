// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test InterpolatedUniformGridData

#include "../tests.h"
#include <deal.II/base/function_lib.h>

// now interpolate the function x*y*z onto points. note that this function is
// (bi/tri)linear and so we can later know what the correct value is that the
// function should provide
Table<1,double> fill (const std_cxx11::array<std::vector<double>,1> &coordinates)
{
  Table<1,double> data(coordinates[0].size());
  for (unsigned int i=0; i<coordinates[0].size(); ++i)
    data[i] = coordinates[0][i];
  return data;
}

Table<2,double> fill (const std_cxx11::array<std::vector<double>,2> &coordinates)
{
  Table<2,double> data(coordinates[0].size(),
		coordinates[1].size());
  for (unsigned int i=0; i<coordinates[0].size(); ++i)
    for (unsigned int j=0; j<coordinates[1].size(); ++j)
      data[i][j] = coordinates[0][i] * coordinates[1][j];
  return data;
}

Table<3,double> fill (const std_cxx11::array<std::vector<double>,3> &coordinates)
{
  Table<3,double> data(coordinates[0].size(),
		coordinates[1].size(),
		coordinates[2].size());
  for (unsigned int i=0; i<coordinates[0].size(); ++i)
    for (unsigned int j=0; j<coordinates[1].size(); ++j)
      for (unsigned int k=0; k<coordinates[2].size(); ++k)
	data[i][j][k] = coordinates[0][i] *
			coordinates[1][j] *
			coordinates[2][k];
  return data;
}


template <int dim>
void check ()
{
  // have coordinate arrays that span an interval starting at d+1
  // d+5 nonuniform intervals
  std_cxx11::array<std::pair<double,double>,dim> intervals;
  std_cxx11::array<unsigned int,dim> n_subintervals;
  for (unsigned int d=0; d<dim; ++d)
    {
      intervals[d] = std::make_pair(d+2., 2*d+5.);
      n_subintervals[d] = d+1 + d*d;
    }

  std_cxx11::array<std::vector<double>,dim> coordinates;
  for (unsigned int d=0; d<dim; ++d)
    {
      const double x = intervals[d].first;
      const double dx = (intervals[d].second-intervals[d].first)/n_subintervals[d];
      
      for (unsigned int i=0; i<n_subintervals[d]+1; ++i)
	coordinates[d].push_back (x+dx*i);
    }
  
  const Table<dim,double> data = fill(coordinates);

  Functions::InterpolatedUniformGridData<dim> f(intervals, n_subintervals, data);

  // now choose a number of randomly chosen points inside the box and
  // verify that the functions returned are correct
  for (unsigned int i=0; i<10; ++i)
    {
      Point<dim> p;
      for (unsigned int d=0; d<dim; ++d)
	p[d] = coordinates[d][0] +
	       (1. * Testing::rand() / RAND_MAX) * (coordinates[d].back() -
						   coordinates[d][0]);

      double exact_value = 1;
      for (unsigned int d=0; d<dim; ++d)
	exact_value *= p[d];
      
      Assert (std::fabs (exact_value - f.value(p)) < 1e-12,
	      ExcInternalError());
    }

  // now also verify that it computes values outside the box correctly, as
  // documented
  double value_at_bottom_left = 1;
  for (unsigned int d=0; d<dim; ++d)
    value_at_bottom_left *= coordinates[d][0];
  
  Assert (std::fabs(f.value(Point<dim>()) - value_at_bottom_left) < 1e-12,
	  ExcInternalError());

  Point<dim> top_right;
  double value_at_top_right = 1;
  for (unsigned int d=0; d<dim; ++d)
    {
      top_right[d] = 1000;
      value_at_top_right *= coordinates[d].back();
    }
  Assert (std::fabs(f.value(top_right) - value_at_top_right) < 1e-12,
	  ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main()
{
  std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1>();
  check<2>();
  check<3>();
}

