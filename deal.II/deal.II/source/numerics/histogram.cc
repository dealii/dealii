//----------------------------  histogram.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  histogram.cc  ---------------------------


#include <base/memory_consumption.h>
#include <lac/vector.h>
#include <numerics/histogram.h>
#include <algorithm>
#include <cmath>


template <typename number>
inline
bool Histogram::logarithmic_less (const number n1,
				  const number n2)
{
  return (((n1<n2) && (n1>0)) ||
	  ((n1<n2) && (n2<=0)) ||
	  ((n2<n1) && (n1>0) && (n2<=0)));
};



Histogram::Interval::Interval (const double left_point,
			       const double right_point) :
		left_point (left_point),
		right_point (right_point),
		content (0)
{};



unsigned int
Histogram::Interval::memory_consumption () const
{
  return sizeof(*this);
};



template <typename number>
void Histogram::evaluate (const typename std::vector<Vector<number> > &values,
			  const std::vector<double>                   &_y_values,
			  const unsigned int                           n_intervals,
			  const IntervalSpacing                        interval_spacing)
{
  Assert (values.size() > 0, ExcEmptyData());
  Assert (n_intervals > 0, ExcInvalidIntervals());
  for (unsigned int i=0; i<values.size(); ++i)
    Assert (values[i].size() > 0, ExcEmptyData());
  Assert (values.size() == _y_values.size(),
	  ExcIncompatibleArraySize(values.size(), _y_values.size()));

				   // store y_values
  y_values = _y_values;
  
				   // first find minimum and maximum value
				   // in the indicators
  number min_value=0, max_value=0;
  switch (interval_spacing)
    {
      case linear:
      {
	min_value = *std::min_element(values[0].begin(),
				      values[0].end());
	max_value = *std::max_element(values[0].begin(),
				      values[0].end());
	
	for (unsigned int i=1; i<values.size(); ++i)
	  {
	    min_value = std::min (min_value,
				  *std::min_element(values[i].begin(),
						    values[i].end()));
	    max_value = std::max (max_value,
				  *std::max_element(values[i].begin(),
						    values[i].end()));
	  };
	
	break;
      };
       
      case logarithmic:
      {
	typedef bool (*comparator) (const number, const number);
	const comparator logarithmic_less_function
	  =
#ifdef __GNU_CC
#if (__GNUC__==2) && (__GNUC_MINOR__ < 95)
	  &logarithmic_less<number>
#else
	  &Histogram::template logarithmic_less<number>
#endif
#else
	  &Histogram::template logarithmic_less<number>	  
#endif
	  ;
	
	min_value = *std::min_element(values[0].begin(),
				      values[0].end(),
				      logarithmic_less_function);
	
	max_value = *std::max_element(values[0].begin(),
				      values[0].end(),
				      logarithmic_less_function);
	
	for (unsigned int i=1; i<values.size(); ++i)
	  {
	    min_value = std::min (min_value,
				  *std::min_element(values[i].begin(),
						    values[i].end(),
						    logarithmic_less_function),
				  logarithmic_less_function);
	    
	    max_value = std::max (max_value,
				  *std::max_element(values[i].begin(),
						    values[i].end(),
						    logarithmic_less_function),
				  logarithmic_less_function);
	  };
	
	break;
      };
       
      default:
	    Assert (false, ExcInternalError());
    };

				   // move right bound arbitrarily if
				   // necessary. sometimes in logarithmic
				   // mode, max_value may be larger than
				   // min_value, but only up to rounding
				   // precision.
  if (max_value <= min_value)
    max_value = min_value+1;


				   // now set up the intervals based on
				   // the min and max values
  intervals.clear ();
				   // set up one list of intervals
				   // for the first data vector. we will
				   // then produce all the other lists
				   // for the other data vectors by
				   // copying
  intervals.push_back (std::vector<Interval>());
  
  switch (interval_spacing)
    {
      case linear:
      {
	const float delta = (max_value-min_value)/n_intervals;

	for (unsigned int n=0; n<n_intervals; ++n)
	  intervals[0].push_back (Interval(min_value+n*delta,
					   min_value+(n+1)*delta));

	break;
      };

      case logarithmic:
      {
	const float delta = (std::log(max_value)-std::log(min_value))/n_intervals;

	for (unsigned int n=0; n<n_intervals; ++n)
	  intervals[0].push_back (Interval(std::exp(std::log(min_value)+n*delta),
					   std::exp(std::log(min_value)+(n+1)*delta)));
	
	break;
      };

      default:
	    Assert (false, ExcInternalError());
    };

				   // fill the other lists of intervals
  for (unsigned int i=1; i<values.size(); ++i)
    intervals.push_back (intervals[0]);


				   // finally fill the intervals
  for (unsigned int i=0; i<values.size(); ++i)
    for (Vector<number>::const_iterator p=values[i].begin();
	 p < values[i].end(); ++p)
      {
					 // find the right place for *p in
					 // intervals[i]. use regular
					 // operator< here instead of
					 // the logarithmic one to
					 // map negative or zero value
					 // to the leftmost interval always
	for (unsigned int n=0; n<n_intervals; ++n)
	  if (*p <= intervals[i][n].right_point)
	    {
	      ++intervals[i][n].content;
	      break;
	    };
      };
};



template <typename number>
void Histogram::evaluate (const Vector<number>    &values,
			  const unsigned int       n_intervals,
			  const IntervalSpacing    interval_spacing) 
{
  std::vector<Vector<number> > values_list (1,
				       values);
  evaluate (values_list, std::vector<double>(1,0.), n_intervals, interval_spacing);
};



void Histogram::write_gnuplot (std::ostream &out) const
{
  AssertThrow (out, ExcIO());
  Assert (intervals.size() > 0, ExcEmptyData());

				   // do a simple 2d plot, if only
				   // one data set is available
  if (intervals.size()==1)
    {
      for (unsigned int n=0; n<intervals[0].size(); ++n)
	out << intervals[0][n].left_point
	    << ' '
	    << intervals[0][n].content
	    << std::endl
	    << intervals[0][n].right_point
	    << ' '
	    << intervals[0][n].content
	    << std::endl;
    }
  else
				     // otherwise create a whole 3d plot
				     // for the data. use th patch method
				     // of gnuplot for this
				     //
				     // run this loop backwards since otherwise
				     // gnuplot thinks the upper side is the
				     // lower side and draws the diagram in
				     // strange colors
    for (int i=intervals.size()-1; i>=0; --i)
      {
	for (unsigned int n=0; n<intervals[i].size(); ++n)
	  out << intervals[i][n].left_point
	      << ' '
	      << (i<static_cast<int>(intervals.size())-1 ?
		  y_values[i+1] :
		  y_values[i] + (y_values[i]-y_values[i-1]))
	      << ' '
	      << intervals[i][n].content
	      << std::endl
	      << intervals[i][n].right_point
	      << ' '
	      << (i<static_cast<int>(intervals.size())-1 ?
		  y_values[i+1] :
		  y_values[i] + (y_values[i]-y_values[i-1]))
	      << ' '
	      << intervals[i][n].content
	      << std::endl;

	out << std::endl;
	for (unsigned int n=0; n<intervals[i].size(); ++n)
	  out << intervals[i][n].left_point
	      << ' '
	      << y_values[i]
	      << ' '
	      << intervals[i][n].content
	      << std::endl
	      << intervals[i][n].right_point
	      << ' '
	      << y_values[i]
	      << ' '
	      << intervals[i][n].content
	      << std::endl;
	
	out << std::endl;
	
      };

  AssertThrow (out, ExcIO());
};



std::string Histogram::get_interval_spacing_names () 
{
  return "linear|logarithmic";
};



Histogram::IntervalSpacing
Histogram::parse_interval_spacing (const std::string &name)
{
  if (name=="linear")
    return linear;
  else
    if (name=="logarithmic")
      return logarithmic;
    else
      {
	AssertThrow (false, ExcInvalidName(name));

	return linear;
      };
};



unsigned int
Histogram::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (intervals) +
	  MemoryConsumption::memory_consumption (y_values));
};



// explicit instantiations for float
template
void Histogram::evaluate (const std::vector<Vector<float> > &values,
			  const std::vector<double>                  &y_values, 
			  const unsigned int                          n_intervals,
			  const IntervalSpacing                       interval_spacing);
template
void Histogram::evaluate (const Vector<float>   &values,
			  const unsigned int     n_intervals,
			  const IntervalSpacing  interval_spacing);


// explicit instantiations for double
template
void Histogram::evaluate (const std::vector<Vector<double> >  &values,
			  const std::vector<double>                    &y_values, 
			  const unsigned int                            n_intervals,
			  const IntervalSpacing                         interval_spacing);
template
void Histogram::evaluate (const Vector<double>   &values,
			  const unsigned int      n_intervals,
			  const IntervalSpacing   interval_spacing);

